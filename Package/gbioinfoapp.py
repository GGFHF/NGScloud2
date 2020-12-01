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
This file contains the classes related to BioInfo application forms in gui mode.
'''

#-------------------------------------------------------------------------------

import os
import PIL.Image
import PIL.ImageTk
import re
import sys
import threading
import tkinter
import tkinter.ttk

import gdialogs
import xbioinfoapp
import xbowtie2
import xbusco
import xcdhit
import xcufflinks
import xconfiguration
import xcutadapt
import xdatabase
import xddradseqtools
import xdetonate
import xec2
import xexpress
import xfastqc
import xgmap
import xhisat2
import xhtseq
import xipyrad
import xkallisto
import xlib
import xngshelper
import xquast
import xraddesigner
import xread
import xreference
import xresult
import xrnaquast
import xsoapdenovo2
import xsoapdenovotrans
import xssh
import xstar
import xstarcode
import xtoa
import xtophat
import xtransabyss
import xtransrate
import xtrimmomatic
import xtrinity

#-------------------------------------------------------------------------------

class FormInstallBioinfoApp(tkinter.Frame):

    #---------------

    def __init__(self, main, app):
        '''
        Execute actions correspending to the creation of a "FormInstallBioinfoApp" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container
        self.app_code = app

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # set the software name
        if self.app_code == xlib.get_bcftools_code():
            self.app_name = xlib.get_bcftools_name()

        elif self.app_code == xlib.get_bedtools_code():
            self.app_name = xlib.get_bedtools_name()

        elif self.app_code == xlib.get_blastplus_code():
            self.app_name = xlib.get_blastplus_name()

        elif self.app_code == xlib.get_bowtie2_code():
            self.app_name = xlib.get_bowtie2_name()

        elif self.app_code == xlib.get_busco_code():
            self.app_name = xlib.get_busco_name()

        elif self.app_code == xlib.get_cd_hit_code():
            self.app_name = xlib.get_cd_hit_name()

        elif self.app_code == xlib.get_cufflinks_code():
            self.app_name = xlib.get_cufflinks_name()

        elif self.app_code == xlib.get_cutadapt_code():
            self.app_name = xlib.get_cutadapt_name()

        elif self.app_code == xlib.get_ddradseqtools_code():
            self.app_name = xlib.get_ddradseqtools_name()

        elif self.app_code == xlib.get_detonate_code():
            self.app_name = xlib.get_detonate_name()

        elif self.app_code == xlib.get_diamond_code():
            self.app_name = xlib.get_diamond_name()

        elif self.app_code == xlib.get_emboss_code():
            self.app_name = xlib.get_emboss_name()

        elif self.app_code == xlib.get_entrez_direct_code():
            self.app_name = xlib.get_entrez_direct_name()

        elif self.app_code == xlib.get_express_code():
            self.app_name = xlib.get_express_name()

        elif self.app_code == xlib.get_fastqc_code():
            self.app_name = xlib.get_fastqc_name()

        elif self.app_code == xlib.get_gmap_gsnap_code():
            self.app_name = xlib.get_gmap_gsnap_name()

        elif self.app_code == xlib.get_hisat2_code():
            self.app_name = xlib.get_hisat2_name()

        elif self.app_code == xlib.get_htseq_code():
            self.app_name = xlib.get_htseq_name()

        elif self.app_code == xlib.get_ipyrad_code():
            self.app_name = xlib.get_ipyrad_name()

        elif self.app_code == xlib.get_kallisto_code():
            self.app_name = xlib.get_kallisto_name()

        elif self.app_code == xlib.get_miniconda3_code():
            self.app_name = xlib.get_miniconda3_name()

        elif self.app_code == xlib.get_ngshelper_code():
            self.app_name = xlib.get_ngshelper_name()

        elif self.app_code == xlib.get_quast_code():
            self.app_name = xlib.get_quast_name()

        elif self.app_code == xlib.get_r_code():
            self.app_name = xlib.get_r_name()

        elif self.app_code == xlib.get_raddesigner_code():
            self.app_name = xlib.get_raddesigner_name()

        elif self.app_code == xlib.get_rnaquast_code():
            self.app_name = xlib.get_rnaquast_name()

        elif self.app_code == xlib.get_rsem_code():
            self.app_name = xlib.get_rsem_name()

        elif self.app_code == xlib.get_samtools_code():
            self.app_name = xlib.get_samtools_name()

        elif self.app_code == xlib.get_soapdenovo2_code():
            self.app_name = xlib.get_soapdenovo2_name()

        elif self.app_code == xlib.get_soapdenovotrans_code():
            self.app_name = xlib.get_soapdenovotrans_name()

        elif self.app_code == xlib.get_star_code():
            self.app_name = xlib.get_star_name()

        elif self.app_code == xlib.get_starcode_code():
            self.app_name = xlib.get_starcode_name()

        elif self.app_code == xlib.get_toa_code():
            self.app_name = xlib.get_toa_name()

        elif self.app_code == xlib.get_tophat_code():
            self.app_name = xlib.get_tophat_name()

        elif self.app_code == xlib.get_transabyss_code():
            self.app_name = xlib.get_transabyss_name()

        elif self.app_code == xlib.get_transdecoder_code():
            self.app_name = xlib.get_transdecoder_name()

        elif self.app_code == xlib.get_transrate_code():
            self.app_name = xlib.get_transrate_name()

        elif self.app_code == xlib.get_trimmomatic_code():
            self.app_name = xlib.get_trimmomatic_name()

        elif self.app_code == xlib.get_trinity_code():
            self.app_name = xlib.get_trinity_name()

        elif self.app_code == xlib.get_vcftools_code():
            self.app_name = xlib.get_vcftools_name()

        elif self.app_code == xlib.get_vcftools_perl_libraries_code():
            self.app_name = xlib.get_vcftools_perl_libraries_name()

        elif self.app_code == xlib.get_vsearch_code():
            self.app_name = xlib.get_vsearch_name()

        # get the version and download URL of the BioInfo application
        (self.bioinfoapp_version, self.bioinfoapp__url, self.bioinfoapp_channels) = xconfiguration.get_bioinfo_app_data(self.app_name)

        # assign the text of the "head"
        self.head = f'{self.app_name} - Install software'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_version_type = tkinter.StringVar()
        self.wrapper_version_type.trace('w', self.check_inputs)
        self.wrapper_version_id = tkinter.StringVar()
        self.wrapper_version_id.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormInstallBioinfoApp".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, columnspan=2, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_version" and register it with the grid geometry manager
        self.label_version = tkinter.Label(self, text='Version')
        self.label_version.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_version_type" and register it with the grid geometry manager
        self.combobox_version_type = tkinter.ttk.Combobox(self, width=13, height=4, textvariable=self.wrapper_version_type)
        self.combobox_version_type.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "entry_version_id" and register it with the grid geometry manager
        self.entry_version_id = tkinter.Entry(self, textvariable=self.wrapper_version_id, width=25, validatecommand=self.check_inputs)
        self.entry_version_id.grid(row=1, column=2, padx=(0,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*66)
        self.label_fit.grid(row=2, column=3, padx=(0,0), pady=(25,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=2, column=5, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_version_type.bind('<<ComboboxSelected>>', self.combobox_version_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.wrapper_version_id.set('')
        if self.bioinfoapp_version in ['any','last']:
            self.wrapper_version_id.set('')
        else:
            self.wrapper_version_id.set(self.bioinfoapp_version)
        self.entry_version_id['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_version_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_version_type(self):
        '''
        Populate data in "combobox_version_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_version_type.set('')

        # load the names of clusters which are running in the combobox
        if self.bioinfoapp_version == 'any':
            self.combobox_version_type['values'] = ['last', 'specific']
            self.combobox_version_type.set('last')
            self.combobox_version_type['state'] = 'readonly'
        elif self.bioinfoapp_version == 'last':
            self.combobox_version_type['values'] = ['last']
            self.combobox_version_type.set('last')
            self.combobox_version_type['state'] = 'disabled'
        else:
            self.combobox_version_type['values'] = ['specific']
            self.combobox_version_type.set('specific')
            self.combobox_version_type['state'] = 'disabled'

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def combobox_version_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_version_type" has been selected
        '''

        # enable or disable the version identification entry
        if self.bioinfoapp_version == 'any':
            if self.wrapper_version_type.get() == 'last':
                self.wrapper_version_id.set('')
                self.entry_version_id['state'] = 'disabled'
            elif self.wrapper_version_type.get() == 'specific':
                self.wrapper_version_id.set('')
                self.entry_version_id['state']='normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormInstallBioinfoApp" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_version_type.get() != '' and (self.wrapper_version_id.get() != '' or self.wrapper_version_type.get() != 'specific'):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the app installation process.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the specific version
        if OK:
            if self.bioinfoapp_version == 'any' and self.wrapper_version_type.get() == 'specific':
                message = f'{xlib.get_project_name()} is designed to use the last version of {self.app_name} and its parameters. The use of some parameter in some previous versions could cause a malfunction..\n\nAre you sure to continue?'
                OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the software installation
        if OK:
            if self.app_code == xlib.get_miniconda3_code():
                message = f'{self.app_name} (Bioconda infrastructure) is going to be installed in the cluster {self.wrapper_cluster_name.get()}. All Bioconda packages previously installed will be lost and they have to be reinstalled.\n\nAre you sure to continue?'
            elif self.app_code == xlib.get_r_code():
                message = f'{self.app_name} and analysis packages are going to be installed in the cluster {self.wrapper_cluster_name.get()}. The previous version will be lost, if it exists.\n\nAre you sure to continue?'
            elif self.app_code in [xlib.get_ddradseqtools_code(), xlib.get_ngshelper_code(), xlib.get_raddesigner_code(), xlib.get_toa_code(), xlib.get_transrate_code()]:
                message = f'{self.app_name} software is going to be installed in the cluster {self.wrapper_cluster_name.get()}. The previous version will be lost, if it exists.\n\nAre you sure to continue?'
            else:
                message = f'The {self.app_name} (channel {self.bioinfoapp_channels}) is going to be installed in the cluster {self.wrapper_cluster_name.get()}. The previous version will be lost, if it exists.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # install the software
        if OK:

            if self.wrapper_version_type.get() == 'last':
                version = 'last'
            else:
                version = self.wrapper_version_id.get()

            # install the BCFTools software
            if self.app_code == xlib.get_bcftools_code():
                package_list = [(xlib.get_bcftools_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the BEDTools software
            elif self.app_code == xlib.get_bedtools_code():
                package_list = [(xlib.get_bedtools_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the BLAST+ software
            elif self.app_code == xlib.get_blastplus_code():
                package_list = [(xlib.get_blastplus_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Bowtie2 software
            elif self.app_code == xlib.get_bowtie2_code():
                package_list = [(xlib.get_bowtie2_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the BUSCO software
            elif self.app_code == xlib.get_busco_code():
                package_list = [(xlib.get_busco_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the CD-HIT software
            elif self.app_code == xlib.get_cd_hit_code():
                package_list = [(xlib.get_cd_hit_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Cufflinks software
            elif self.app_code == xlib.get_cufflinks_code():
                package_list = [(xlib.get_cufflinks_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the cutadapt software
            elif self.app_code == xlib.get_cutadapt_code():
                package_list = [(xlib.get_cutadapt_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the ddRADseqTools software
            elif self.app_code == xlib.get_ddradseqtools_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.install_ddradseqtools.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.install_ddradseqtools, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the DETONATE software
            elif self.app_code == xlib.get_detonate_code():
                package_list = [(xlib.get_detonate_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_bowtie2_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the DIAMOND software
            elif self.app_code == xlib.get_diamond_code():
                package_list = [(xlib.get_diamond_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the EMBOSS software
            elif self.app_code == xlib.get_emboss_code():
                package_list = [(xlib.get_emboss_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Entrez Direct software
            elif self.app_code == xlib.get_entrez_direct_code():
                package_list = [(xlib.get_entrez_direct_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the eXpress software
            elif self.app_code == xlib.get_express_code():
                package_list = [(xlib.get_express_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the FastQC software
            elif self.app_code == xlib.get_fastqc_code():
                package_list = [(xlib.get_fastqc_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the GMAP-GSNAP software
            elif self.app_code == xlib.get_gmap_gsnap_code():
                package_list = [(xlib.get_gmap_gsnap_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the HISAT2 software
            elif self.app_code == xlib.get_hisat2_code():
                package_list = [(xlib.get_hisat2_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the HTSeq software
            elif self.app_code == xlib.get_htseq_code():
                package_list = [(xlib.get_htseq_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the ipyrad software
            elif self.app_code == xlib.get_ipyrad_code():
                package_list = [(xlib.get_ipyrad_conda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the kallisto software
            elif self.app_code == xlib.get_kallisto_code():
                package_list = [(xlib.get_kallisto_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Miniconda3 software
            elif self.app_code == xlib.get_miniconda3_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_miniconda3.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_miniconda3, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the NGShelper software
            elif self.app_code == xlib.get_ngshelper_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xngshelper.install_ngshelper.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xngshelper.install_ngshelper, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the QUAST software
            elif self.app_code == xlib.get_quast_code():
                package_list = [(xlib.get_quast_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install R and analysis packages
            elif self.app_code == xlib.get_r_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_r.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_r, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the RADdesigner software
            elif self.app_code == xlib.get_raddesigner_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xraddesigner.install_raddesigner.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xraddesigner.install_raddesigner, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the rnaQUAST software
            elif self.app_code == xlib.get_rnaquast_code():
                package_list = [(xlib.get_rnaquast_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the RSEM software
            elif self.app_code == xlib.get_rsem_code():
                package_list = [(xlib.get_rsem_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the SAMtools software
            if self.app_code == xlib.get_samtools_code():
                package_list = [(xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the SOAPdenovo2 software
            elif self.app_code == xlib.get_soapdenovo2_code():
                package_list = [(xlib.get_soapdenovo2_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the SOAPdenovo-Trans software
            elif self.app_code == xlib.get_soapdenovotrans_code():
                package_list = [(xlib.get_soapdenovotrans_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the STAR software
            elif self.app_code == xlib.get_star_code():
                package_list = [(xlib.get_star_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the starcode software
            elif self.app_code == xlib.get_starcode_code():
                package_list = [(xlib.get_starcode_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the TOA software
            elif self.app_code == xlib.get_toa_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtoa.install_toa.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtoa.install_toa, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the TopHat software
            elif self.app_code == xlib.get_tophat_code():
                package_list = [(xlib.get_tophat_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_samtools_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_tabix_anaconda_code(), 'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Trans-ABySS software
            elif self.app_code == xlib.get_transabyss_code():
                package_list = [(xlib.get_transabyss_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the TransDecoder software
            elif self.app_code == xlib.get_transdecoder_code():
                package_list = [(xlib.get_transdecoder_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Transrate software
            elif self.app_code == xlib.get_transrate_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtransrate.install_transrate.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtransrate.install_transrate, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()
                # -- package_list = [(xlib.get_transrate_anaconda_code(), version, self.bioinfoapp_channels)]
                # -- dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                # -- threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                # -- threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Trimmomatic software
            elif self.app_code == xlib.get_trimmomatic_code():
                package_list = [(xlib.get_trimmomatic_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the Trinity software
            elif self.app_code == xlib.get_trinity_code():
                package_list = [(xlib.get_trinity_anaconda_code(), version, self.bioinfoapp_channels), (xlib.get_bowtie2_anaconda_code() ,'last', self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the VCFtools software
            elif self.app_code == xlib.get_vcftools_code():
                package_list = [(xlib.get_vcftools_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the VCFtools Perl libraries software
            elif self.app_code == xlib.get_vcftools_perl_libraries_code():
                package_list = [(xlib.get_vcftools_perl_libraries_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # install the vsearch software
            elif self.app_code == xlib.get_vsearch_code():
                package_list = [(xlib.get_vsearch_anaconda_code(), version, self.bioinfoapp_channels)]
                dialog_log = gdialogs.DialogLog(self, self.head, xbioinfoapp.install_anaconda_package_list.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbioinfoapp.install_anaconda_package_list, args=(self.app_code, self.app_name, package_list, self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormInstallBioinfoApp".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateBowtie2ConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateBowtie2ConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_bowtie2_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateBowtie2ConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=7, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('NONE')
        self.combobox_assembly_dataset['state'] = 'disabled'
        self.assembly_dataset_id = 'NONE'
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('NONE')
        self.combobox_assembly_type['state'] = 'disabled'
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code(), xlib.get_soapdenovo2_code(), xlib.get_starcode_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file", "combobox_assembly_dataset" and "combobox_assembly_type"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
            self.assembly_dataset_id = ''
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = []
            self.wrapper_assembly_dataset.set('')
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = []
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'disabled'

        else:
            self.populate_combobox_reference_file()
            self.assembly_dataset_id = 'NONE'
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = ['NONE']
            self.wrapper_assembly_dataset.set('NONE')
            self.combobox_assembly_dataset['state'] = 'disabled'
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.populate_combobox_assembly_dataset()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or self.assembly_dataset_id.startswith(xlib.get_soapdenovo2_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()) or self.assembly_dataset_id.startswith(xlib.get_starcode_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateBowtie2ConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the Bowtie2 config file
        if OK:
            message = f'The file {xbowtie2.get_bowtie2_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Bowtie2 config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xbowtie2.create_bowtie2_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xbowtie2.create_bowtie2_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Bowtie2 config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xbowtie2.get_bowtie2_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xbowtie2.check_bowtie2_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_bowtie2_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateBowtie2ConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateBuscoConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateBuscoConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_busco_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateBuscoConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*49)
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateBuscoConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != ''  and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the BUSCO config file
        if OK:
            message = f'The file {xbusco.get_busco_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the BUSCO config file
        if OK:
            (OK, error_list) = xbusco.create_busco_config_file(self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the BUSCO config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xbusco.get_busco_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xbusco.check_busco_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_busco_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateBuscoConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateCdHitEstConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateCdHitEstConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_cd_hit_est_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateCdHitEstConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*49)
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateCdHitEstConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the CD-HIT-EST config file
        if OK:
            message = f'The file {xcdhit.get_cd_hit_est_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the CD-HIT-EST config file
        if OK:
            (OK, error_list) = xcdhit.create_cd_hit_est_config_file(self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the CD-HIT-EST config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xcdhit.get_cd_hit_est_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xcdhit.check_cd_hit_est_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_cd_hit_est_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateCdHitEstConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateCuffdiffCuffnormConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main, app):
        '''
        Execute actions correspending to the creation of a "FormRecreateCuffdiffCuffnormConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container
        self.app_code = app

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # set the software name
        if self.app_code == xlib.get_cuffdiff_code():
            self.app_name = xlib.get_cuffdiff_name()
        elif self.app_code == xlib.get_cuffnorm_code():
            self.app_name = xlib.get_cuffnorm_name()

        # assign the text of the "head"
        self.head = f'{self.app_name} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_quantitation_dataset = tkinter.StringVar()
        self.wrapper_quantitation_dataset.trace('w', self.check_inputs)
        self.wrapper_abundance_files = tkinter.StringVar()
        self.wrapper_abundance_files.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateCuffdiffCuffnormConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_select_dirs"
        image_select_dirs = PIL.Image.open('./image_select_dirs.png')
        imagetk_select_dirs = PIL.ImageTk.PhotoImage(image_select_dirs)

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_quantitation_dataset" and register it with the grid geometry manager
        self.label_quantitation_dataset = tkinter.Label(self, text='Quantitation dataset')
        self.label_quantitation_dataset.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_quantitation_dataset" and register it with the grid geometry manager
        self.combobox_quantitation_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_quantitation_dataset)
        self.combobox_quantitation_dataset.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_abundance_files" and register it with the grid geometry manager
        self.label_abundance_files = tkinter.Label(self, text='Abundance files')
        self.label_abundance_files.grid(row=4, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "entry_abundance_files" and register it with the grid geometry manager
        self.entry_abundance_files = tkinter.Entry(self, textvariable=self.wrapper_abundance_files, width=45, state='disabled', validatecommand=self.check_inputs)
        self.entry_abundance_files.grid(row=4, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "button_select_abundance_files" and register it with the grid geometry manager
        self.button_select_abundance_files = tkinter.ttk.Button(self, image=imagetk_select_dirs, command=self.select_abundance_files, state='disabled')
        self.button_select_abundance_files.image = imagetk_select_dirs
        self.button_select_abundance_files.grid(row=4, column=2, padx=(5,0), pady=(40,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*46)
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_quantitation_dataset.bind('<<ComboboxSelected>>', self.combobox_quantitation_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_quantitation_dataset['values'] = []
        self.wrapper_quantitation_dataset.set('')
        self.quantitation_dataset_id = None
        self.wrapper_abundance_files.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_cufflinks_cuffmerge_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def populate_combobox_quantitation_dataset(self):
        '''
        Populate data in "combobox_quantitation_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_quantitation_dataset.set('')

        # get the list of the quantitation dataset names
        app_list = [xlib.get_cuffquant_code()]
        (_, _, quantitation_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the quantitation dataset names in the combobox
        self.combobox_quantitation_dataset['values'] = sorted(quantitation_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_assembly_dataset()

        # load data in "combobox_quantitation_dataset"
        self.populate_combobox_quantitation_dataset()

        # initialize the data of alignment datasets
        self.wrapper_abundance_files.set('')
        self.abundance_file_id_list = []
        self.group_list = []

        # enable "button_select_abundance_files"
        self.button_select_abundance_files['state'] = 'disabled'

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_quantitation_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_quantitation_dataset" has been selected
        '''

        # get the quantitation dataset identification
        (_, _, self.quantitation_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_quantitation_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # get the path of quantitation dataset in the cluster
        cluster_quantitation_dataset_id = xlib.get_cluster_experiment_result_dataset_dir(self.wrapper_experiment_id.get(), self.quantitation_dataset_id)

        # initialize the data of alignment datasets
        self.wrapper_abundance_files.set('')
        self.abundance_file_id_list = []
        self.group_list = []

        # enable "button_select_abundance_files"
        self.button_select_abundance_files['state'] = 'enabled'

        # get the abundance file list
        command = f'cd {cluster_quantitation_dataset_id}; find . -type f -regex "./.*cxb"'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                self.abundance_file_id_list.append(os.path.basename(line.rstrip('\n')))
        else:
            message = f'*** ERROR: Wrong command ---> {command}'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
        if self.abundance_file_id_list == []:
            message = f'WARNING: There are not abundance files in the cluster directory {self.quantitation_dataset_id}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
        else:
            self.abundance_file_id_list.sort()
            for i in range(len(self.abundance_file_id_list)):
                self.group_list.append(f'group-{i+1}')

    #---------------

    def select_abundance_files(self):
        '''
        Set the group of the abundance file and update "entry_abundance_files".
        '''

        # get the directory dictionary of directories in the volume
        abundance_file_dict = {}
        for i in range(len(self.abundance_file_id_list)):
            key = self.abundance_file_id_list[i]
            abundance_file_dict[key] = {'option_id': self.abundance_file_id_list[i], 'option_value': self.group_list[i], 'comment': '', 'value_type': 'any', 'admitted_option_value_list': []}

        # build the data dictionary
        data_dict = {}
        data_dict['option_id']= {'text': 'Library abundance file', 'width': 60, 'alignment': 'left'}
        data_dict['option_value'] = {'text': 'Group', 'width': 10, 'alignment': 'left'}
        data_dict['comment'] = {'text': 'Admitted values', 'width': 0, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        title_text = 'Update the group to which the libraries belong'
        window_height = 600
        window_width = 575
        auxliary_window_height = 60
        auxliary_window_width = 920
        dialog_table = gdialogs.DialogOptionUpdate(self, title_text, window_height, window_width, auxliary_window_height, auxliary_window_width, data_dict, abundance_file_dict, self.abundance_file_id_list)
        self.wait_window(dialog_table)

        # update the group list
        for key in abundance_file_dict.keys():
            index = self.abundance_file_id_list.index(key)
            self.group_list[index] = abundance_file_dict[key]['option_value']

        # update "entry_abundance_files"
        self.wrapper_abundance_files.set(str(self.abundance_file_id_list).strip('[]').replace('\'',''))

    #---------------

    def check_inputs(self, *args):
        '''
        check the content of each input of "FormRecreateCuffdiffCuffnormConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_quantitation_dataset.get() != '' and self.wrapper_abundance_files.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Cuffdiff config file
        if OK:
            if self.app_code == xlib.get_cuffdiff_code():
                message = f'The file {xcufflinks.get_cuffdiff_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            elif self.app_code == xlib.get_cuffnorm_code():
                message = f'The file {xcufflinks.get_cuffnorm_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Cuffdiff/Cuffnorm config file
        if OK:
            if self.app_code == xlib.get_cuffdiff_code():
                (OK, error_list) = xcufflinks.create_cuffdiff_config_file(self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.quantitation_dataset_id, self.abundance_file_id_list, self.group_list)
            elif self.app_code == xlib.get_cuffnorm_code():
                (OK, error_list) = xcufflinks.create_cuffnorm_config_file(self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.quantitation_dataset_id, self.abundance_file_id_list, self.group_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Cuffdiff/Cuffnorm config file
        if OK:

            # edit the config file using "DialogEditor" 
            if self.app_code == xlib.get_cuffdiff_code():
                dialog_editor = gdialogs.DialogEditor(self, xcufflinks.get_cuffdiff_config_file())
            elif self.app_code == xlib.get_cuffnorm_code():
                dialog_editor = gdialogs.DialogEditor(self, xcufflinks.get_cuffnorm_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            if self.app_code == xlib.get_cuffdiff_code():
                (OK, error_list) = xcufflinks.check_cuffdiff_config_file(strict=False)
            elif self.app_code == xlib.get_cuffnorm_code():
                (OK, error_list) = xcufflinks.check_cuffnorm_config_file(strict=False)
            if OK:
                message = f'The {self.app_name} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateCuffdiffCuffnormConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateCufflinksCuffmergeConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateCufflinksCuffmergeConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_cufflinks_cuffmerge_name()} - Recreate config file'

        # initialize the selected alignment dataset list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_mask_file = tkinter.StringVar()
        self.wrapper_mask_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_alignment_datasets = tkinter.StringVar()
        self.wrapper_alignment_datasets.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateCufflinksCuffmergeConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_select_dirs"
        image_select_dirs = PIL.Image.open('./image_select_dirs.png')
        imagetk_select_dirs = PIL.ImageTk.PhotoImage(image_select_dirs)

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_mask_file" and register it with the grid geometry manager
        self.label_mask_file = tkinter.Label(self, text='Mask file')
        self.label_mask_file.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_mask_file" and register it with the grid geometry manager
        self.combobox_mask_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_mask_file)
        self.combobox_mask_file.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_alignment_datasets" and register it with the grid geometry manager
        self.label_alignment_datasets = tkinter.Label(self, text='Alignment datasets')
        self.label_alignment_datasets.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_alignment_datasets" and register it with the grid geometry manager
        self.entry_alignment_datasets = tkinter.Entry(self, textvariable=self.wrapper_alignment_datasets, width=45, state='disabled', validatecommand=self.check_inputs)
        self.entry_alignment_datasets.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "button_select_alignment_datasets" and register it with the grid geometry manager
        self.button_select_alignment_datasets = tkinter.ttk.Button(self, image=imagetk_select_dirs, command=self.select_alignment_datasets, state='disabled')
        self.button_select_alignment_datasets.image = imagetk_select_dirs
        self.button_select_alignment_datasets.grid(row=6, column=2, padx=(5,0), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*36)
        self.label_fit.grid(row=7, column=3, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=4, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=5, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_mask_file.bind('<<ComboboxSelected>>', self.combobox_mask_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_mask_file['values'] = []
        self.wrapper_mask_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.wrapper_alignment_datasets.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the annotation file names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the annotation file names in the combobox
        self.combobox_annotation_file['values'] = sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_mask_file(self):
        '''
        Populate data in "combobox_mask_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_mask_file.set('')

        # get the list of the mask file names
        (_, _, mask_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the mask file names in the combobox
        self.combobox_mask_file['values'] = ['NONE'] + sorted(mask_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)


    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # clear data in "combobox_mask_file"
        self.combobox_mask_file['values'] = []
        self.wrapper_mask_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # load data in "combobox_annotation_file"
        self.populate_combobox_annotation_file()

        # load data in "combobox_mask_file"
        self.populate_combobox_mask_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_mask_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_mask_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # initialize the data of alignment datasets
        self.wrapper_alignment_datasets.set('')
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # enable "button_select_alignment_datasets"
        self.button_select_alignment_datasets['state'] = 'enable'

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def select_alignment_datasets(self):
        '''
        Select the alignment dataset identifications and update "entry_alignment_datasets".
        '''

        # get the list of the alignment dataset names
        app_list = xcufflinks.get_alignment_software_code_list()
        (OK, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # get the directory dictionary of directories in the volume
        if OK:
            alignment_dataset_dict = {}
            for alignment_dataset_name in alignment_dataset_name_list:
                key = alignment_dataset_name
                option_value = 'YES' if alignment_dataset_name in self.selected_alignment_dataset_list else 'NO'
                alignment_dataset_dict[key] = {'option_id': alignment_dataset_name, 'option_value': option_value, 'comment': 'YES or NO', 'value_type': 'uppercase_string_list', 'admitted_option_value_list': ['YES', 'NO']}

        # build the data dictionary
        if OK:
            data_dict = {}
            data_dict['option_id']= {'text': 'Alignment dataset', 'width': 22, 'alignment': 'left'}
            data_dict['option_value'] = {'text': 'Is selected?', 'width': 12, 'alignment': 'left'}
            data_dict['comment'] = {'text': 'Admitted values', 'width': 17, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        if OK:
            title_text = 'Select alignment datasets'
            window_height = 600
            window_width = 420
            auxliary_window_height = 60
            auxliary_window_width = 760
            dialog_table = gdialogs.DialogOptionUpdate(self, title_text, window_height, window_width, auxliary_window_height, auxliary_window_width, data_dict, alignment_dataset_dict, sorted(alignment_dataset_dict.keys()))
            self.wait_window(dialog_table)

        # recreate the selected alignment dataset list and alignment dataset identification list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []
        for key in alignment_dataset_dict.keys():
            if alignment_dataset_dict[key]['option_value'] == 'YES':
                self.selected_alignment_dataset_list.append(key)
                (OK, _, alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), key, status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)
                self.alignment_dataset_id_list.append(alignment_dataset_id)

        # sort the selected alignment dataset list and alignment dataset identification list
        if self.selected_alignment_dataset_list != []:
            self.selected_alignment_dataset_list.sort()
            self.alignment_dataset_id_list.sort()

        # update "entry_alignment_datasets"
        self.wrapper_alignment_datasets.set(str(self.selected_alignment_dataset_list).strip('[]').replace('\'',''))

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateCufflinksCuffmergeConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_annotation_file.get() != '' and self.wrapper_mask_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.alignment_dataset_id_list != []:
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Cufflinks-Cuffmerge config file
        if OK:
            message = f'The file {xcufflinks.get_cufflinks_cuffmerge_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Cufflinks-Cuffmerge config file
        if OK:
            (OK, error_list) = xcufflinks.create_cufflinks_cuffmerge_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.wrapper_mask_file.get(), self.alignment_dataset_id_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Cufflinks-Cuffmerge config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xcufflinks.get_cufflinks_cuffmerge_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xcufflinks.check_cufflinks_cuffmerge_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_cufflinks_cuffmerge_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateCufflinksCuffmergeConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateCuffquantConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateCuffquantConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_cuffquant_name()} - Recreate config file'

        # initialize the selected alignment dataset list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_mask_file = tkinter.StringVar()
        self.wrapper_mask_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_alignment_datasets = tkinter.StringVar()
        self.wrapper_alignment_datasets.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateCuffquantConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_select_dirs"
        image_select_dirs = PIL.Image.open('./image_select_dirs.png')
        imagetk_select_dirs = PIL.ImageTk.PhotoImage(image_select_dirs)

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_mask_file" and register it with the grid geometry manager
        self.label_mask_file = tkinter.Label(self, text='Mask file')
        self.label_mask_file.grid(row=2, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_mask_file" and register it with the grid geometry manager
        self.combobox_mask_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_mask_file)
        self.combobox_mask_file.grid(row=2, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_alignment_datasets" and register it with the grid geometry manager
        self.label_alignment_datasets = tkinter.Label(self, text='Alignment datasets')
        self.label_alignment_datasets.grid(row=4, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "entry_alignment_datasets" and register it with the grid geometry manager
        self.entry_alignment_datasets = tkinter.Entry(self, textvariable=self.wrapper_alignment_datasets, width=45, state='disabled', validatecommand=self.check_inputs)
        self.entry_alignment_datasets.grid(row=4, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "button_select_alignment_datasets" and register it with the grid geometry manager
        self.button_select_alignment_datasets = tkinter.ttk.Button(self, image=imagetk_select_dirs, command=self.select_alignment_datasets, state='disabled')
        self.button_select_alignment_datasets.image = imagetk_select_dirs
        self.button_select_alignment_datasets.grid(row=4, column=2, padx=(5,0), pady=(40,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=5, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=5, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*36)
        self.label_fit.grid(row=6, column=3, padx=(0,0), pady=(40,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=6, column=4, padx=(0,5), pady=(40,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=6, column=5, padx=(5,5), pady=(40,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_mask_file.bind('<<ComboboxSelected>>', self.combobox_mask_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_mask_file['values'] = []
        self.wrapper_mask_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.wrapper_alignment_datasets.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_mask_file(self):
        '''
        Populate data in "combobox_mask_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_mask_file.set('')

        # get the list of the mask file names
        (_, _, mask_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the mask file names in the combobox
        self.combobox_mask_file['values'] = ['NONE'] + sorted(mask_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_cufflinks_cuffmerge_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_mask_file"
        self.combobox_mask_file['values'] = []
        self.wrapper_mask_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_mask_file"
        self.populate_combobox_mask_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_mask_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_mask_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # initialize the data of alignment datasets
        self.wrapper_alignment_datasets.set('')
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # enable "button_select_alignment_datasets"
        self.button_select_alignment_datasets['state'] = 'enable'

        # load data in "combobox_read_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def select_alignment_datasets(self):
        '''
        Select the alignment dataset identifications and update "entry_alignment_datasets".
        '''

        # get the list of the alignment dataset names
        app_list = xcufflinks.get_alignment_software_code_list()

        (OK, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # get the directory dictionary of directories in the volume
        if OK:
            alignment_dataset_dict = {}
            for alignment_dataset_name in alignment_dataset_name_list:
                key = alignment_dataset_name
                option_value = 'YES' if alignment_dataset_name in self.selected_alignment_dataset_list else 'NO'
                alignment_dataset_dict[key] = {'option_id': alignment_dataset_name, 'option_value': option_value, 'comment': 'YES or NO', 'value_type': 'uppercase_string_list', 'admitted_option_value_list': ['YES', 'NO']}

        # build the data dictionary
        if OK:
            data_dict = {}
            data_dict['option_id']= {'text': 'Alignment dataset', 'width': 22, 'alignment': 'left'}
            data_dict['option_value'] = {'text': 'Is selected?', 'width': 12, 'alignment': 'left'}
            data_dict['comment'] = {'text': 'Admitted values', 'width': 17, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        if OK:
            title_text = 'Select alignment datasets'
            window_height = 600
            window_width = 420
            auxliary_window_height = 60
            auxliary_window_width = 760
            dialog_table = gdialogs.DialogOptionUpdate(self, title_text, window_height, window_width, auxliary_window_height, auxliary_window_width, data_dict, alignment_dataset_dict, sorted(alignment_dataset_dict.keys()))
            self.wait_window(dialog_table)

        # recreate the selected alignment dataset list and alignment dataset identification list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []
        for key in alignment_dataset_dict.keys():
            if alignment_dataset_dict[key]['option_value'] == 'YES':
                self.selected_alignment_dataset_list.append(key)
                (OK, _, alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), key, status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)
                self.alignment_dataset_id_list.append(alignment_dataset_id)

        # sort the selected alignment dataset list and alignment dataset identification list
        if self.selected_alignment_dataset_list != []:
            self.selected_alignment_dataset_list.sort()
            self.alignment_dataset_id_list.sort()

        # update "entry_alignment_datasets"
        self.wrapper_alignment_datasets.set(str(self.selected_alignment_dataset_list).strip('[]').replace('\'',''))

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateCuffquantConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_mask_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.alignment_dataset_id_list != [] and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Cuffquant config file
        if OK:
            message = f'The file {xcufflinks.get_cuffquant_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Cuffquant config file
        if OK:
            (OK, error_list) = xcufflinks.create_cuffquant_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(),  self.wrapper_mask_file.get(), self.alignment_dataset_id_list, self.assembly_dataset_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Cuffquant config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xcufflinks.get_cuffquant_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xcufflinks.check_cuffquant_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_cuffquant_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateCuffquantConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateCutadaptConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateCutadaptConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_cutadapt_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateCutadaptConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateCutadaptConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the cutadapt config file
        if OK:
            message = f'The file {xcutadapt.get_cutadapt_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the cutadapt config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xcutadapt.create_cutadapt_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xcutadapt.create_cutadapt_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the cutadapt config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xcutadapt.get_cutadapt_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xcutadapt.check_cutadapt_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_cutadapt_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateCutadaptConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateDdRadSeqSimulationConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateDdRadSeqSimulationConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_ddradseq_simulation_name()} - Recreate config file'

        # initialize the enzyme identification list
        self.enzyme_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_enzyme1 = tkinter.StringVar()
        self.wrapper_enzyme1.trace('w', self.check_inputs)
        self.wrapper_enzyme2 = tkinter.StringVar()
        self.wrapper_enzyme2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateDdRadSeqSimulationConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme1" and register it with the grid geometry manager
        self.label_enzyme1 = tkinter.Label(self, text='Enzyme1')
        self.label_enzyme1.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_enzyme1" and register it with the grid geometry manager
        self.combobox_enzyme1 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme1)
        self.combobox_enzyme1.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme1_warning" and register it with the grid geometry manager
        self.label_enzyme1_warning = tkinter.Label(self, text='')
        self.label_enzyme1_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme2" and register it with the grid geometry manager
        self.label_enzyme2 = tkinter.Label(self, text='Enzyme2')
        self.label_enzyme2.grid(row=4, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_enzyme2" and register it with the grid geometry manager
        self.combobox_enzyme2 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme2)
        self.combobox_enzyme2.grid(row=4, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme2_warning" and register it with the grid geometry manager
        self.label_enzyme2_warning = tkinter.Label(self, text='')
        self.label_enzyme2_warning.grid(row=4, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*49)
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_enzyme1.bind('<<ComboboxSelected>>', self.combobox_enzyme1_selected_item)
        self.combobox_enzyme1.bind('<FocusOut>', self.combobox_enzyme1_focus_out)
        self.combobox_enzyme2.bind('<<ComboboxSelected>>', self.combobox_enzyme2_selected_item)
        self.combobox_enzyme2.bind('<FocusOut>', self.combobox_enzyme2_focus_out)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_enzyme1['values'] = []
        self.wrapper_enzyme1.set('')
        self.enzyme1_id = None
        self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)
        self.combobox_enzyme2['values'] = []
        self.wrapper_enzyme2.set('')
        self.enzyme2_id = None
        self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_enzyme1(self):
        '''
        Populate data in "combobox_enzyme1".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme1.set('')

        # load the ids of enzymes
        self.combobox_enzyme1['values'] = self.enzyme_id_seq_list

    #---------------

    def populate_combobox_enzyme2(self):
        '''
        Populate data in "combobox_enzyme2".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme2.set('')

        # load the ids of enzymes
        self.combobox_enzyme2['values'] = self.enzyme_id_seq_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # get the dictionary of restriction enzymes
        (OK, error_list, self.restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        if not OK:
            message = ''
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            self.close()

        # get the list of enzyme identification and restriction site seq
        if OK:
            self.enzyme_id_seq_list = xddradseqtools.get_enzyme_id_seq_list(self.restriction_enzyme_dict)

        # get the enzime identification list
        if OK:
            self.enzyme_id_list = list(self.restriction_enzyme_dict.keys())
            self.enzyme_id_list.sort()

        # load data in "combobox_enzyme1"
        if OK:
            self.combobox_enzyme1['state'] = 'normal'
            self.populate_combobox_enzyme1()

        # load data in "combobox_enzyme2"
        if OK:
            self.combobox_enzyme2['state'] = 'normal'
            self.populate_combobox_enzyme2()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_enzyme1_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme1" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

    #---------------

    def combobox_enzyme1_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme1" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

        # check that "combobox_enzyme1" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme1_id, self.label_enzyme1_warning)

    #---------------

    def combobox_enzyme2_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme2" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

    #---------------

    def combobox_enzyme2_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme2" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

        # check that "combobox_enzyme2" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme2_id, self.label_enzyme2_warning)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateDdRadSeqSimulationConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_enzyme1.get() != '' and self.wrapper_enzyme2.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_enzyme(self, enzyme, enzyme_warning):
        '''
        Check the content enzyme
        '''

        # initialize the control variable
        OK = True

        # check that enzyme value is an identifier of a restriction enzyme or a restriction site sequence
        if enzyme != '' and enzyme not in self.enzyme_id_seq_list and enzyme not in self.enzyme_id_list and not xlib.is_valid_sequence(seq=enzyme, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
            enzyme_warning['text'] = 'Invalid enzyme id or restriction site seq.'
            enzyme_warning['foreground'] = 'red'
            OK = False
        else:
            enzyme_warning['text'] = 'Enzyme id or restriction site seq (cut point: *).'
            enzyme_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs() and self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)  and self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check that enzyme 1 has to be different from enzyme 2
        if OK:
            if self.combobox_enzyme1.get() in self.enzyme_id_list:
                enzyme1_seq = self.restriction_enzyme_dict[self.combobox_enzyme1.get()]['restriction_site_seq']
            else:
                id1 = xddradseqtools.get_enzyme_id(self.combobox_enzyme1.get(), self.restriction_enzyme_dict)
                if id1 != None:
                    enzyme1_seq = self.restriction_enzyme_dict[id1]['restriction_site_seq']
                else:
                    enzyme1_seq = self.combobox_enzyme1.get()
            if self.combobox_enzyme2.get() in self.enzyme_id_list:
                enzyme2_seq = self.restriction_enzyme_dict[self.combobox_enzyme2.get()]['restriction_site_seq']
            else:
                id2 = xddradseqtools.get_enzyme_id(self.combobox_enzyme2.get(), self.restriction_enzyme_dict)
                if id2 != None:
                    enzyme2_seq = self.restriction_enzyme_dict[id2]['restriction_site_seq']
                else:
                    enzyme2_seq = self.combobox_enzyme2.get()
            if enzyme1_seq.upper() == enzyme2_seq.upper():
                OK = False
                message = 'Both enzymes have the same sequence. A ddRADseq experiment has to be performed with two different enzymes.'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the ddRADseq simulation config file
        if OK:
            message = f'The file {xddradseqtools.get_ddradseq_simulation_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the ddRADseq simulation config file
        if OK:
            (OK, error_list) = xddradseqtools.create_ddradseq_simulation_config_file(self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.enzyme1_id, self.enzyme2_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the ddRADseq simulation config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xddradseqtools.get_ddradseq_simulation_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xddradseqtools.check_ddradseq_simulation_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_ddradseq_simulation_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateDdRadSeqSimulationConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateExpressConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateExpressConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_express_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # initialize the selected alignment dataset list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)
        self.wrapper_alignment_datasets = tkinter.StringVar()
        self.wrapper_alignment_datasets.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateExpressConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_select_dirs"
        image_select_dirs = PIL.Image.open('./image_select_dirs.png')
        imagetk_select_dirs = PIL.ImageTk.PhotoImage(image_select_dirs)

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_alignment_datasets" and register it with the grid geometry manager
        self.label_alignment_datasets = tkinter.Label(self, text='Alignment datasets')
        self.label_alignment_datasets.grid(row=4, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_alignment_datasets" and register it with the grid geometry manager
        self.entry_alignment_datasets = tkinter.Entry(self, textvariable=self.wrapper_alignment_datasets, width=45, state='disabled', validatecommand=self.check_inputs)
        self.entry_alignment_datasets.grid(row=4, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "button_select_alignment_datasets" and register it with the grid geometry manager
        self.button_select_alignment_datasets = tkinter.ttk.Button(self, image=imagetk_select_dirs, command=self.select_alignment_datasets, state='disabled')
        self.button_select_alignment_datasets.image = imagetk_select_dirs
        self.button_select_alignment_datasets.grid(row=4, column=2, padx=(5,0), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*50)
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'
        self.wrapper_alignment_datasets.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # initialize the data of alignment datasets
        self.wrapper_alignment_datasets.set('')
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # enable "button_select_alignment_datasets"
        self.button_select_alignment_datasets['state'] = 'enable'

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def select_alignment_datasets(self):
        '''
        Select the alignment dataset identifications and update "entry_alignment_datasets".
        '''

        # get the list of the alignment dataset names
        app_list = xexpress.get_alignment_software_code_list()
        (OK, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # get the directory dictionary of directories in the volume
        if OK:
            alignment_dataset_dict = {}
            for alignment_dataset_name in alignment_dataset_name_list:
                key = alignment_dataset_name
                option_value = 'YES' if alignment_dataset_name in self.selected_alignment_dataset_list else 'NO'
                alignment_dataset_dict[key] = {'option_id': alignment_dataset_name, 'option_value': option_value, 'comment': 'YES or NO', 'value_type': 'uppercase_string_list', 'admitted_option_value_list': ['YES', 'NO']}

        # build the data dictionary
        if OK:
            data_dict = {}
            data_dict['option_id']= {'text': 'Alignment dataset', 'width': 22, 'alignment': 'left'}
            data_dict['option_value'] = {'text': 'Is selected?', 'width': 12, 'alignment': 'left'}
            data_dict['comment'] = {'text': 'Admitted values', 'width': 17, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        if OK:
            title_text = 'Select alignment datasets'
            window_height = 600
            window_width = 420
            auxliary_window_height = 60
            auxliary_window_width = 760
            dialog_table = gdialogs.DialogOptionUpdate(self, title_text, window_height, window_width, auxliary_window_height, auxliary_window_width, data_dict, alignment_dataset_dict, sorted(alignment_dataset_dict.keys()))
            self.wait_window(dialog_table)

        # recreate the selected alignment dataset list and alignment dataset identification list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []
        for key in alignment_dataset_dict.keys():
            if alignment_dataset_dict[key]['option_value'] == 'YES':
                self.selected_alignment_dataset_list.append(key)
                (_, _, alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), key, status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)
                self.alignment_dataset_id_list.append(alignment_dataset_id)

        # sort the selected alignment dataset list and alignment dataset identification list
        if self.selected_alignment_dataset_list != []:
            self.selected_alignment_dataset_list.sort()
            self.alignment_dataset_id_list.sort()

        # update "entry_alignment_datasets"
        self.wrapper_alignment_datasets.set(str(self.selected_alignment_dataset_list).strip('[]').replace('\'',''))

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateExpressConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != ''  and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type != '' and self.alignment_dataset_id_list != []:
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Express config file
        if OK:
            message = f'The file {xexpress.get_express_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Express config file
        if OK:
            (OK, error_list) = xexpress.create_express_config_file(self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.alignment_dataset_id_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Express config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xexpress.get_express_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xexpress.check_express_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_express_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateExpressConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateFastQCConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateFastQCConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_fastqc_name()} - Recreate config file'

        # initialize the read dataset identification
        self.read_dataset_id = None

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateFastQCConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*55)
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.wrapper_file_pattern.set('.*fastq')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_run_set"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        check the content of each input of "FormRecreateFastQCConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False

        # confirm the creation of the FastQC config file
        if OK:
            message = f'The file {xfastqc.get_fastqc_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the FastQC config file
        if OK:
            (OK, error_list) = xfastqc.create_fastqc_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, selected_file_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the FastQC config file corresponding to the environment
        if OK:
            
            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xfastqc.get_fastqc_config_file())
            self.wait_window(dialog_editor)
            
            # check the config file
            (OK, error_list) = xfastqc.check_fastqc_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_fastqc_name()} transfer config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateFastQCConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateGgTrinityConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateGgTrinityConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_ggtrinity_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_alignment_dataset = tkinter.StringVar()
        self.wrapper_alignment_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateGgTrinityConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_alignment_dataset" and register it with the grid geometry manager
        self.label_alignment_dataset = tkinter.Label(self, text='Alignment dataset')
        self.label_alignment_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_alignment_dataset" and register it with the grid geometry manager
        self.combobox_alignment_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_alignment_dataset)
        self.combobox_alignment_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*48)
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_alignment_dataset.bind('<<ComboboxSelected>>', self.combobox_alignment_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_alignment_dataset['values'] = []
        self.wrapper_alignment_dataset.set('')
        self.alignment_dataset_id = None

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_alignment_dataset(self):
        '''
        Populate data in "combobox_alignment_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_alignment_dataset.set('')

        # get the list of the alignment dataset names
        app_list = [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code(), xlib.get_star_code(), xlib.get_tophat_code()]
        (_, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the alignment dataset names in the combobox
        self.combobox_alignment_dataset['values'] = sorted(alignment_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_alignment_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_alignment_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_alignment_dataset" has been selected
        '''

        # get the alignment dataset identification
        (_, _, self.alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_alignment_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateGgTrinityConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_alignment_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Genome-guided Trinity config file
        if OK:
            message = f'The file {xtrinity.get_ggtrinity_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Cuffdiff config file
        if OK:
            (OK, error_list) = xtrinity.create_ggtrinity_config_file(self.wrapper_experiment_id.get(), self.alignment_dataset_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Genome-guided Trinity config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtrinity.get_ggtrinity_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtrinity.check_ggtrinity_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_ggtrinity_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateGgTrinityConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateGmapConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateGmapConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_gmap_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateGmapConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*48)
        self.label_fit.grid(row=6, column=2, padx=(0,0), pady=(40,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=6, column=3, padx=(0,5), pady=(40,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=6, column=4, padx=(5,5), pady=(40,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateGmapConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the GMAP config file
        if OK:
            message = f'The file {xgmap.get_gmap_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the GMAP config file
        if OK:
            (OK, error_list) = xgmap.create_gmap_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the GMAP config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xgmap.get_gmap_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xgmap.check_gmap_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_gmap_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateGmapConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateGsnapConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateGsnapConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_gsnap_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateGsnapConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=7, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('NONE')
        self.combobox_assembly_dataset['state'] = 'disabled'
        self.assembly_dataset_id = 'NONE'
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('NONE')
        self.combobox_assembly_type['state'] = 'disabled'
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code(), xlib.get_soapdenovo2_code(), xlib.get_starcode_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file", "combobox_assembly_dataset" and "combobox_assembly_type"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
            self.assembly_dataset_id = ''
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = []
            self.wrapper_assembly_dataset.set('')
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = []
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'disabled'

        else:
            self.populate_combobox_reference_file()
            self.assembly_dataset_id = 'NONE'
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = ['NONE']
            self.wrapper_assembly_dataset.set('NONE')
            self.combobox_assembly_dataset['state'] = 'disabled'
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.populate_combobox_assembly_dataset()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or self.assembly_dataset_id.startswith(xlib.get_soapdenovo2_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()) or self.assembly_dataset_id.startswith(xlib.get_starcode_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateGsnapConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the GSNAP config file
        if OK:
            message = f'The file {xgmap.get_gsnap_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the GSNAP config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xgmap.create_gsnap_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xgmap.create_gsnap_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the GSNAP config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xgmap.get_gsnap_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xgmap.check_gsnap_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_gsnap_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateGsnapConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateHisat2ConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateHisat2ConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_hisat2_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateHisat2ConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=6, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=10, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=10, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=10, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the annotation file names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the annotation file names in the combobox
        self.combobox_annotation_file['values'] = ['NONE'] + sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()


        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # load data in "combobox_annotation_file"
        self.populate_combobox_annotation_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateHisat2ConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != ''  and self.wrapper_annotation_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the HISAT2 config file
        if OK:
            message = f'The file {xhisat2.get_hisat2_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the HISAT2 config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xhisat2.create_hisat2_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xhisat2.create_hisat2_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the HISAT2 config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xhisat2.get_hisat2_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xhisat2.check_hisat2_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_hisat2_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateHisat2ConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateHtseqCountConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateHtseqCountConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_htseq_count_name()} - Recreate config file'

        # initialize the selected alignment dataset list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_alignment_datasets = tkinter.StringVar()
        self.wrapper_alignment_datasets.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateHtseqCountConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_select_dirs"
        image_select_dirs = PIL.Image.open('./image_select_dirs.png')
        imagetk_select_dirs = PIL.ImageTk.PhotoImage(image_select_dirs)

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_alignment_datasets" and register it with the grid geometry manager
        self.label_alignment_datasets = tkinter.Label(self, text='Alignment datasets')
        self.label_alignment_datasets.grid(row=4, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_alignment_datasets" and register it with the grid geometry manager
        self.entry_alignment_datasets = tkinter.Entry(self, textvariable=self.wrapper_alignment_datasets, width=45, state='disabled', validatecommand=self.check_inputs)
        self.entry_alignment_datasets.grid(row=4, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "button_select_alignment_datasets" and register it with the grid geometry manager
        self.button_select_alignment_datasets = tkinter.ttk.Button(self, image=imagetk_select_dirs, command=self.select_alignment_datasets, state='disabled')
        self.button_select_alignment_datasets.image = imagetk_select_dirs
        self.button_select_alignment_datasets.grid(row=4, column=2, padx=(5,0), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*48)
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.wrapper_alignment_datasets.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the annotation file names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the annotation file names in the combobox
        self.combobox_annotation_file['values'] = sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)


    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_annotation_file"
        self.populate_combobox_annotation_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # initialize the data of alignment datasets
        self.wrapper_alignment_datasets.set('')
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []

        # enable "button_select_alignment_datasets"
        self.button_select_alignment_datasets['state'] = 'enable'

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def select_alignment_datasets(self):
        '''
        Select the alignment dataset identifications and update "entry_alignment_datasets".
        '''

        # get the list of the alignment dataset names
        app_list = xhtseq.get_alignment_software_code_list()
        (OK, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # get the directory dictionary of directories in the volume
        if OK:
            alignment_dataset_dict = {}
            for alignment_dataset_name in alignment_dataset_name_list:
                key = alignment_dataset_name
                option_value = 'YES' if alignment_dataset_name in self.selected_alignment_dataset_list else 'NO'
                alignment_dataset_dict[key] = {'option_id': alignment_dataset_name, 'option_value': option_value, 'comment': 'YES or NO', 'value_type': 'uppercase_string_list', 'admitted_option_value_list': ['YES', 'NO']}

        # build the data dictionary
        if OK:
            data_dict = {}
            data_dict['option_id']= {'text': 'Alignment dataset', 'width': 22, 'alignment': 'left'}
            data_dict['option_value'] = {'text': 'Is selected?', 'width': 12, 'alignment': 'left'}
            data_dict['comment'] = {'text': 'Admitted values', 'width': 17, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        if OK:
            title_text = 'Select alignment datasets'
            window_height = 600
            window_width = 420
            auxliary_window_height = 60
            auxliary_window_width = 760
            dialog_table = gdialogs.DialogOptionUpdate(self, title_text, window_height, window_width, auxliary_window_height, auxliary_window_width, data_dict, alignment_dataset_dict, sorted(alignment_dataset_dict.keys()))
            self.wait_window(dialog_table)

        # recreate the selected alignment dataset list and alignment dataset identification list
        self.selected_alignment_dataset_list = []
        self.alignment_dataset_id_list = []
        for key in alignment_dataset_dict.keys():
            if alignment_dataset_dict[key]['option_value'] == 'YES':
                self.selected_alignment_dataset_list.append(key)
                (_, _, alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), key, status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)
                self.alignment_dataset_id_list.append(alignment_dataset_id)

        # sort the selected alignment dataset list and alignment dataset identification list
        if self.selected_alignment_dataset_list != []:
            self.selected_alignment_dataset_list.sort()
            self.alignment_dataset_id_list.sort()

        # update "entry_alignment_datasets"
        self.wrapper_alignment_datasets.set(str(self.selected_alignment_dataset_list).strip('[]').replace('\'',''))

    #---------------

    def check_inputs(self, *args):
        '''
        check the content of each input of "FormRecreateHtseqCountConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_annotation_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.alignment_dataset_id_list != []:
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the htseq-count config file
        if OK:
            message = f'The file {xhtseq.get_htseq_count_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the htseq-count config file
        if OK:
            (OK, error_list) = xhtseq.create_htseq_count_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_annotation_file.get(), self.alignment_dataset_id_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the htseq-count config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xhtseq.get_htseq_count_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xhtseq.check_htseq_count_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_htseq_count_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateHtseqCountConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateInsilicoReadNormalizationConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateInsilicoReadNormalizationConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_insilico_read_normalization_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateInsilicoReadNormalizationConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(5,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateInsilicoReadNormalizationConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the insilico_read_normalization config file
        if OK:
            message = f'The file {xtrinity.get_insilico_read_normalization_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the insilico_read_normalization config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtrinity.create_insilico_read_normalization_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xtrinity.create_insilico_read_normalization_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the insilico_read_normalization config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtrinity.get_insilico_read_normalization_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtrinity.check_insilico_read_normalization_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_insilico_read_normalization_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateInsilicoReadNormalizationConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateIpyradConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateIpyradConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_ipyrad_name()} - Recreate config file'

        # initialize the enzyme identification list
        self.enzyme_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_assembly_method = tkinter.StringVar()
        self.wrapper_assembly_method.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_datatype = tkinter.StringVar()
        self.wrapper_datatype.trace('w', self.check_inputs)
        self.wrapper_enzyme1 = tkinter.StringVar()
        self.wrapper_enzyme1.trace('w', self.check_inputs)
        self.wrapper_enzyme2 = tkinter.StringVar()
        self.wrapper_enzyme2.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateIpyradConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_assembly_method" and register it with the grid geometry manager
        self.label_assembly_method = tkinter.Label(self, text='Assembly method')
        self.label_assembly_method.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_method" and register it with the grid geometry manager
        self.combobox_assembly_method = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_assembly_method)
        self.combobox_assembly_method.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='disable', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='disable', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_datatype" and register it with the grid geometry manager
        self.label_datatype = tkinter.Label(self, text='RAD datatype')
        self.label_datatype.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_datatype" and register it with the grid geometry manager
        self.combobox_datatype = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_datatype)
        self.combobox_datatype.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_enzyme1" and register it with the grid geometry manager
        self.label_enzyme1 = tkinter.Label(self, text='Enzyme1')
        self.label_enzyme1.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_enzyme1" and register it with the grid geometry manager
        self.combobox_enzyme1 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme1)
        self.combobox_enzyme1.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_enzyme1_warning" and register it with the grid geometry manager
        self.label_enzyme1_warning = tkinter.Label(self, text='')
        self.label_enzyme1_warning.grid(row=5, column=2, columnspan= 3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_enzyme2" and register it with the grid geometry manager
        self.label_enzyme2 = tkinter.Label(self, text='Enzyme2')
        self.label_enzyme2.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_enzyme2" and register it with the grid geometry manager
        self.combobox_enzyme2 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme2)
        self.combobox_enzyme2.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_enzyme2_warning" and register it with the grid geometry manager
        self.label_enzyme2_warning = tkinter.Label(self, text='')
        self.label_enzyme2_warning.grid(row=6, column=2, columnspan= 3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=9, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*50)
        self.label_fit.grid(row=10, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=10, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=10, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_assembly_method.bind('<<ComboboxSelected>>', self.combobox_assembly_method_selected_item)
        self.combobox_datatype.bind('<<ComboboxSelected>>', self.combobox_datatype_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_enzyme1.bind('<<ComboboxSelected>>', self.combobox_enzyme1_selected_item)
        self.combobox_enzyme1.bind('<FocusOut>', self.combobox_enzyme1_focus_out)
        self.combobox_enzyme2.bind('<<ComboboxSelected>>', self.combobox_enzyme2_selected_item)
        self.combobox_enzyme2.bind('<FocusOut>', self.combobox_enzyme2_focus_out)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_assembly_method['values'] = []
        self.wrapper_assembly_method.set('')
        self.combobox_datatype['values'] = []
        self.wrapper_datatype.set('')
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_enzyme1['values'] = []
        self.wrapper_enzyme1.set('')
        self.enzyme1_id = None
        self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)
        self.combobox_enzyme2['values'] = []
        self.wrapper_enzyme2.set('')
        self.enzyme2_id = None
        self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_assembly_method(self):
        '''
        Populate data in "combobox_assembly_method".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_method.set('')

        # load the list of the assembly methods in the combobox
        self.combobox_assembly_method['values'] =['DENOVO', 'REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE']

    #---------------

    def populate_combobox_datatype(self):
        '''
        Populate data in "combobox_datatype".
        '''

        # clear the value selected in the combobox
        self.wrapper_datatype.set('')

        # load the list of the assembly methods in the combobox
        self.combobox_datatype['values'] =['RAD', 'DDRAD', 'PAIRDDRAD', 'GBS', 'PAIRGBS']

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_enzyme1(self):
        '''
        Populate data in "combobox_enzyme1".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme1.set('')

        # load the ids of enzymes
        self.combobox_enzyme1['values'] = self.enzyme_id_seq_list

    #---------------

    def populate_combobox_enzyme2(self):
        '''
        Populate data in "combobox_enzyme2".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme2.set('')

        # load the ids of enzymes
        self.combobox_enzyme2['values'] = self.enzyme_id_seq_list

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # get the dictionary of restriction enzymes
        (OK, error_list, self.restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        if not OK:
            message = ''
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            self.close()

        # load data in "combobox_assembly_method"
        self.populate_combobox_assembly_method()

        # load data in "combobox_datatype"
        self.populate_combobox_datatype()

        # get the list of enzyme identification and restriction site seq
        if OK:
            self.enzyme_id_seq_list = xddradseqtools.get_enzyme_id_seq_list(self.restriction_enzyme_dict)

        # get the enzime identification list
        if OK:
            self.enzyme_id_list = list(self.restriction_enzyme_dict.keys())
            self.enzyme_id_list.sort()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_method_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_method" has been selected
        '''

        # if the method is DENOVO
        if self.wrapper_assembly_method.get() == 'DENOVO':

            # set NONE in "combobox_reference_dataset"
            self.combobox_reference_dataset.set('NONE')
            self.combobox_reference_dataset['state'] = 'disabled'

            # set NONE in "combobox_reference_file"
            self.combobox_reference_file.set('NONE')
            self.combobox_reference_file['state'] = 'disabled'

        # if the method is REFERENCE, DENOVO+REFERENCE or DENOVO-REFERENCE
        elif self.wrapper_assembly_method.get() in ['REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE']:

            # load data in "combobox_reference_dataset"
            self.combobox_reference_dataset['state'] = 'readonly'
            self.populate_combobox_reference_dataset()

            # clear data in "combobox_reference_file"
            self.combobox_reference_file['state'] = 'readonly'
            self.combobox_reference_file['values'] = []
            self.wrapper_reference_file.set('')

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_datatype_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_datatype" has been selected
        '''

        # if the method is DDRAD or PAIRDDRAD
        if self.wrapper_datatype.get() in ['DDRAD', 'PAIRDDRAD']:

            # load data in "combobox_enzyme1"
            self.combobox_enzyme1['state'] = 'normal'
            self.populate_combobox_enzyme1()

            # load data in "combobox_enzyme2"
            self.combobox_enzyme2['state'] = 'normal'
            self.populate_combobox_enzyme2()

        # if the method is other one
        else:

            # load data in "combobox_enzyme1"
            self.combobox_enzyme1['state'] = 'normal'
            self.populate_combobox_enzyme1()

            # set NONE in "combobox_enzyme2"
            self.combobox_enzyme2.set('NONE')
            self.combobox_enzyme2['state'] = 'disabled'

    #---------------

    def combobox_enzyme1_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme1" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

    #---------------

    def combobox_enzyme1_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme1" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

        # check that "combobox_enzyme1" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme1_id, self.label_enzyme1_warning)

    #---------------

    def combobox_enzyme2_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme2" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

    #---------------

    def combobox_enzyme2_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme2" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

        # check that "combobox_enzyme2" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme2_id, self.label_enzyme2_warning)

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateIpyradConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_assembly_method.get() != '' and self.wrapper_datatype.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_enzyme1.get() != '' and self.wrapper_enzyme2.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def check_enzyme(self, enzyme, enzyme_warning):
        '''
        Check the content enzyme
        '''

        # initialize the control variable
        OK = True

        # check that enzyme value is an identifier of a restriction enzyme or a restriction site sequence
        if enzyme != '' and enzyme not in self.enzyme_id_seq_list and enzyme not in self.enzyme_id_list and not xlib.is_valid_sequence(seq=enzyme, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
            enzyme_warning['text'] = 'Invalid enzyme id or restriction site seq.'
            enzyme_warning['foreground'] = 'red'
            OK = False
        else:
            enzyme_warning['text'] = 'Enzyme id or restriction site seq (cut point: *).'
            enzyme_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        if self.wrapper_datatype.get() in ['DDRAD', 'PAIRDDRAD']:
            OK = self.check_inputs() and self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)  and self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)
        else:
            OK = self.check_inputs() and self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check that enzyme 1 has to be different from enzyme 2 when the method is DDRAD or PAIRDDRAD
        if OK:
            if self.wrapper_datatype.get() in ['DDRAD', 'PAIRDDRAD']:
                if self.combobox_enzyme1.get() in self.enzyme_id_list:
                    enzyme1_seq = self.restriction_enzyme_dict[self.combobox_enzyme1.get()]['restriction_site_seq']
                else:
                    id1 = xddradseqtools.get_enzyme_id(self.combobox_enzyme1.get(), self.restriction_enzyme_dict)
                    if id1 != None:
                        enzyme1_seq = self.restriction_enzyme_dict[id1]['restriction_site_seq']
                    else:
                        enzyme1_seq = self.combobox_enzyme1.get()
                if self.combobox_enzyme2.get() in self.enzyme_id_list:
                    enzyme2_seq = self.restriction_enzyme_dict[self.combobox_enzyme2.get()]['restriction_site_seq']
                else:
                    id2 = xddradseqtools.get_enzyme_id(self.combobox_enzyme2.get(), self.restriction_enzyme_dict)
                    if id2 != None:
                        enzyme2_seq = self.restriction_enzyme_dict[id2]['restriction_site_seq']
                    else:
                        enzyme2_seq = self.combobox_enzyme2.get()
                if enzyme1_seq.upper() == enzyme2_seq.upper():
                    OK = False
                    message = 'Both enzymes have the same sequence. A ddRADseq experiment has to be performed with two different enzymes.'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the ipyrad config file
        if OK:
            message = f'The file {xipyrad.get_ipyrad_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the ipyrad config file
        if OK:
            (OK, error_list) = xipyrad.create_ipyrad_config_file(self.wrapper_experiment_id.get(), self.wrapper_assembly_method.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_datatype.get(), self.enzyme1_id, self.enzyme2_id, self.read_dataset_id, self.wrapper_file_pattern.get(), 'NO')
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the ipyrad config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xipyrad.get_ipyrad_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xipyrad.check_ipyrad_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_ipyrad_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateIpyradConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateKallistoConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateKallistoConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_kallisto_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateKallistoConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=5, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the reference dataset names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_annotation_file['values'] = sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_annotation_file"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_annotation_file['values'] = ['NONE']
            self.wrapper_annotation_file.set('NONE')
        else:
            self.populate_combobox_annotation_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateKallistoConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_annotation_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != '') and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the kallisto config file
        if OK:
            message = f'The file {xkallisto.get_kallisto_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the kallisto config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xkallisto.create_kallisto_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            elif self.read_type == 'PE':
                (OK, error_list) = xkallisto.create_kallisto_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the kallisto config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xkallisto.get_kallisto_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xkallisto.check_kallisto_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_kallisto_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateKallistoConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateQuastConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateQuastConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_quast_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateQuastConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*50)
        self.label_fit.grid(row=6, column=2, padx=(0,0), pady=(40,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=6, column=3, padx=(0,5), pady=(40,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=6, column=4, padx=(5,5), pady=(40,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # load data in "combobox_reference_file"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
        else:
            self.populate_combobox_reference_file()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateQuastConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the QUAST config file
        if OK:
            message = f'The file {xquast.get_quast_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the QUAST config file
        if OK:
            (OK, error_list) = xquast.create_quast_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the QUAST config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xquast.get_quast_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xquast.check_quast_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_quast_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateQuastConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateRADdesignerConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateRADdesignerConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_raddesigner_name()} - Recreate config file'

        # initialize the read dataset identification
        self.vcf_location_dataset_id = None

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_vcf_location_dataset = tkinter.StringVar()
        self.wrapper_vcf_location_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateRADdesignerConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_vcf_location_dataset" and register it with the grid geometry manager
        self.label_vcf_location_dataset = tkinter.Label(self, text='VCF location dataset')
        self.label_vcf_location_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_vcf_location_dataset" and register it with the grid geometry manager
        self.combobox_vcf_location_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_vcf_location_dataset)
        self.combobox_vcf_location_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*46)
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_vcf_location_dataset.bind('<<ComboboxSelected>>', self.combobox_vcf_location_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_vcf_location_dataset['values'] = []
        self.wrapper_vcf_location_dataset.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_vcf_location_dataset(self):
        '''
        Populate data in "combobox_vcf_location_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_vcf_location_dataset.set('')

        # get the list of the read dataset names
        app_list = [xlib.get_ipyrad_code()]
        (_, _, vcf_location_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_vcf_location_dataset['values'] = vcf_location_dataset_name_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_vcf_location_dataset"
        self.combobox_vcf_location_dataset['values'] = []
        self.wrapper_vcf_location_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_run_set"
        self.populate_combobox_vcf_location_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_vcf_location_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read dataset identification
        (_, _, self.vcf_location_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_vcf_location_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        check the content of each input of "FormRecreateRADdesignerConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_vcf_location_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the RADdesigner config file
        if OK:
            message = f'The file {xraddesigner.get_raddesigner_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the RADdesigner config file
        if OK:
            (OK, error_list) = xraddesigner.create_raddesigner_config_file(self.wrapper_experiment_id.get(), self.vcf_location_dataset_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the RADdesigner config file corresponding to the environment
        if OK:
            
            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xraddesigner.get_raddesigner_config_file())
            self.wait_window(dialog_editor)
            
            # check the config file
            (OK, error_list) = xraddesigner.check_raddesigner_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_raddesigner_name()} transfer config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateRADdesignerConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateRefEvalConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateRefEvaltConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_ref_eval_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateRefEvaltConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=5, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*18)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
        else:
            self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateRefEvaltConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != '') and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the REF-EVAL config file
        if OK:
            message = f'The file {xdetonate.get_ref_eval_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the REF-EVAL config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xdetonate.create_ref_eval_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            elif self.read_type == 'PE':
                (OK, error_list) = xdetonate.create_ref_eval_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the REF-EVAL config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xdetonate.get_ref_eval_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xdetonate.check_ref_eval_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_ref_eval_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateRefEvaltConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateRnaQuastConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateRnaQuastConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_rnaquast_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateRnaQuastConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=5, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
        else:
            self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateRnaQuastConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != '') and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the rnaQUAST config file
        if OK:
            message = f'The file {xrnaquast.get_rnaquast_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the rnaQUAST config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xrnaquast.create_rnaquast_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            elif self.read_type == 'PE':
                (OK, error_list) = xrnaquast.create_rnaquast_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the rnaQUAST config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xrnaquast.get_rnaquast_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xrnaquast.check_rnaquast_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_rnaquast_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateRnaQuastConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateRsemEvalConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateRsemEvalConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_rsem_eval_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateRsemEvalConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=7, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=7, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=8, column=0, padx=(15,5), pady=(20,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=8, column=1, padx=(5,5), pady=(20,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=9, column=2, padx=(0,0), pady=(20,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=9, column=3, padx=(0,5), pady=(20,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=9, column=4, padx=(5,5), pady=(20,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateRsemEvalConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != '') and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the RSEM-EVAL config file
        if OK:
            message = f'The file {xdetonate.get_rsem_eval_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the RSEM-EVAL config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xdetonate.create_rsem_eval_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            elif self.read_type == 'PE':
                (OK, error_list) = xdetonate.create_rsem_eval_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list, self.assembly_dataset_id, self.wrapper_assembly_type.get())                 
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the RSEM-EVAL config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xdetonate.get_rsem_eval_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xdetonate.check_rsem_eval_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_rsem_eval_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateRsemEvalConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateRSiteSearchConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateRSiteSearchConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_rsitesearch_name()} - Recreate config file'

        # initialize the enzyme identification list
        self.enzyme_id_list = []

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_enzyme1 = tkinter.StringVar()
        self.wrapper_enzyme1.trace('w', self.check_inputs)
        self.wrapper_enzyme2 = tkinter.StringVar()
        self.wrapper_enzyme2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateRSiteSearchConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme1" and register it with the grid geometry manager
        self.label_enzyme1 = tkinter.Label(self, text='Enzyme1')
        self.label_enzyme1.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_enzyme1" and register it with the grid geometry manager
        self.combobox_enzyme1 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme1)
        self.combobox_enzyme1.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme1_warning" and register it with the grid geometry manager
        self.label_enzyme1_warning = tkinter.Label(self, text='')
        self.label_enzyme1_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme2" and register it with the grid geometry manager
        self.label_enzyme2 = tkinter.Label(self, text='Enzyme2')
        self.label_enzyme2.grid(row=4, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_enzyme2" and register it with the grid geometry manager
        self.combobox_enzyme2 = tkinter.ttk.Combobox(self, width=20, height=4, state='disabled', textvariable=self.wrapper_enzyme2)
        self.combobox_enzyme2.grid(row=4, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_enzyme2_warning" and register it with the grid geometry manager
        self.label_enzyme2_warning = tkinter.Label(self, text='')
        self.label_enzyme2_warning.grid(row=4, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*49)
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_enzyme1.bind('<<ComboboxSelected>>', self.combobox_enzyme1_selected_item)
        self.combobox_enzyme1.bind('<FocusOut>', self.combobox_enzyme1_focus_out)
        self.combobox_enzyme2.bind('<<ComboboxSelected>>', self.combobox_enzyme2_selected_item)
        self.combobox_enzyme2.bind('<FocusOut>', self.combobox_enzyme2_focus_out)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_enzyme1['values'] = []
        self.wrapper_enzyme1.set('')
        self.enzyme1_id = None
        self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)
        self.combobox_enzyme2['values'] = []
        self.wrapper_enzyme2.set('')
        self.enzyme2_id = None
        self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_enzyme1(self):
        '''
        Populate data in "combobox_enzyme1".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme1.set('')

        # load the ids of enzymes
        self.combobox_enzyme1['values'] = self.enzyme_id_seq_list

    #---------------

    def populate_combobox_enzyme2(self):
        '''
        Populate data in "combobox_enzyme2".
        '''

        # clear the value selected in the combobox
        self.wrapper_enzyme2.set('')

        # load the ids of enzymes
        self.combobox_enzyme2['values'] = self.enzyme_id_seq_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # get the dictionary of restriction enzymes
        (OK, error_list, self.restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
        if not OK:
            message = ''
            for error in error_list:
                message = f'{message}{error}\n'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            self.close()

        # get the list of enzyme identification and restriction site seq
        if OK:
            self.enzyme_id_seq_list = xddradseqtools.get_enzyme_id_seq_list(self.restriction_enzyme_dict)

        # get the enzime identification list
        if OK:
            self.enzyme_id_list = list(self.restriction_enzyme_dict.keys())
            self.enzyme_id_list.sort()

        # load data in "combobox_enzyme1"
        if OK:
            self.combobox_enzyme1['state'] = 'normal'
            self.populate_combobox_enzyme1()

        # load data in "combobox_enzyme2"
        if OK:
            self.combobox_enzyme2['state'] = 'normal'
            self.populate_combobox_enzyme2()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_enzyme1_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme1" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

    #---------------

    def combobox_enzyme1_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme1" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme1.get() in self.enzyme_id_seq_list:
            self.enzyme1_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme1.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme1_id = self.wrapper_enzyme1.get()

        # check that "combobox_enzyme1" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme1_id, self.label_enzyme1_warning)

    #---------------

    def combobox_enzyme2_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_enzyme2" has been selected
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

    #---------------

    def combobox_enzyme2_focus_out(self, event=None):
        '''
        Process the event when the focus was moved from "combobox_enzyme2" to another widget
        '''

        # get the enzyme identification
        if self.wrapper_enzyme2.get() in self.enzyme_id_seq_list:
            self.enzyme2_id = xddradseqtools.get_enzyme_id(self.wrapper_enzyme2.get(), self.restriction_enzyme_dict)
        else:
            self.enzyme2_id = self.wrapper_enzyme2.get()

        # check that "combobox_enzyme2" value is an identifier of a restriction enzyme or a restriction site sequence
        self.check_enzyme(self.enzyme2_id, self.label_enzyme2_warning)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateRSiteSearchConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_enzyme1.get() != '' and self.wrapper_enzyme2.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_enzyme(self, enzyme, enzyme_warning):
        '''
        Check the content enzyme
        '''

        # initialize the control variable
        OK = True

        # check that enzyme value is an identifier of a restriction enzyme or a restriction site sequence
        if enzyme != '' and enzyme not in self.enzyme_id_seq_list and enzyme not in self.enzyme_id_list and not xlib.is_valid_sequence(seq=enzyme, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
            enzyme_warning['text'] = 'Invalid enzyme id or restriction site seq.'
            enzyme_warning['foreground'] = 'red'
            OK = False
        else:
            enzyme_warning['text'] = 'Enzyme id or restriction site seq (cut point: *).'
            enzyme_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs() and self.check_enzyme(self.combobox_enzyme1.get(), self.label_enzyme1_warning)  and self.check_enzyme(self.combobox_enzyme2.get(), self.label_enzyme2_warning)
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check if enzyme 1 is different or not to enzyme 2
        if OK:
            if self.combobox_enzyme1.get() in self.enzyme_id_list:
                enzyme1_seq = self.restriction_enzyme_dict[self.combobox_enzyme1.get()]['restriction_site_seq']
            else:
                id1 = xddradseqtools.get_enzyme_id(self.combobox_enzyme1.get(), self.restriction_enzyme_dict)
                if id1 != None:
                    enzyme1_seq = self.restriction_enzyme_dict[id1]['restriction_site_seq']
                else:
                    enzyme1_seq = self.combobox_enzyme1.get()
            if self.combobox_enzyme2.get() in self.enzyme_id_list:
                enzyme2_seq = self.restriction_enzyme_dict[self.combobox_enzyme2.get()]['restriction_site_seq']
            else:
                id2 = xddradseqtools.get_enzyme_id(self.combobox_enzyme2.get(), self.restriction_enzyme_dict)
                if id2 != None:
                    enzyme2_seq = self.restriction_enzyme_dict[id2]['restriction_site_seq']
                else:
                    enzyme2_seq = self.combobox_enzyme2.get()
            if enzyme1_seq.upper() == enzyme2_seq.upper():
                message = 'Both enzymes have the same sequence. An enzyme analysis of a RAD-seq experiment will be performed.'
            else:
                message = 'Both enzymes have different sequences. An enzyme analysis of a ddRADseq experiment will be performed.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the rsitesearch config file
        if OK:
            message = f'The file {xddradseqtools.get_rsitesearch_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the rsitesearch config file
        if OK:
            (OK, error_list) = xddradseqtools.create_rsitesearch_config_file(self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.enzyme1_id, self.enzyme2_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the rsitesearch config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xddradseqtools.get_rsitesearch_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xddradseqtools.check_rsitesearch_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_rsitesearch_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateRSiteSearchConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateSoapdenovo2ConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateSoapdenovo2ConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_soapdenovo2_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateSoapdenovo2ConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateSoapdenovo2ConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the SOAPdenovo2 config file
        if OK:
            message = f'The file {xsoapdenovo2.get_soapdenovo2_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the SOAPdenovo2 config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xsoapdenovo2.create_soapdenovo2_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xsoapdenovo2.create_soapdenovo2_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the SOAPdenovo2 config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xsoapdenovo2.get_soapdenovo2_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xsoapdenovo2.check_soapdenovo2_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_soapdenovo2_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateSoapdenovo2ConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateSoapdenovoTransConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateSoapdenovoTransConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_soapdenovotrans_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateSoapdenovoTransConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateSoapdenovoTransConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the SOAPdenovo-Trans config file
        if OK:
            message = f'The file {xsoapdenovotrans.get_soapdenovotrans_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the SOAPdenovo-Trans config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xsoapdenovotrans.create_soapdenovotrans_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xsoapdenovotrans.create_soapdenovotrans_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the SOAPdenovo-Trans config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xsoapdenovotrans.get_soapdenovotrans_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xsoapdenovotrans.check_soapdenovotrans_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_soapdenovotrans_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateSoapdenovoTransConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateSTARConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateSTARConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_star_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateSTARConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=6, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=10, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=10, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=10, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the annotation file names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the annotation file names in the combobox
        self.combobox_annotation_file['values'] = sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # load data in "combobox_annotation_file"
        self.populate_combobox_annotation_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateSTARConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_annotation_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the STAR config file
        if OK:
            message = f'The file {xstar.get_star_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the STAR config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xstar.create_star_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xstar.create_star_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the STAR config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xstar.get_star_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xstar.check_star_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_star_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateSTARConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateStarcodeConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateStarcodeConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_starcode_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateStarcodeConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        # -- self.combobox_read_type['values'] =['Single-end', 'Paired-end']
        self.combobox_read_type['values'] =['Single-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateStarcodeConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the starcode config file
        if OK:
            message = f'The file {xstarcode.get_starcode_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the starcode config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xstarcode.create_starcode_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xstarcode.create_starcode_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the starcode config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xstarcode.get_starcode_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xstarcode.check_starcode_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_starcode_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateStarcodeConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTopHatConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTopHatConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_tophat_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_annotation_file = tkinter.StringVar()
        self.wrapper_annotation_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTopHatConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(50,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(50,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_annotation_file" and register it with the grid geometry manager
        self.label_annotation_file = tkinter.Label(self, text='Annotation file')
        self.label_annotation_file.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_annotation_file" and register it with the grid geometry manager
        self.combobox_annotation_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_annotation_file)
        self.combobox_annotation_file.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=6, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=10, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=10, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=10, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_annotation_file.bind('<<ComboboxSelected>>', self.combobox_annotation_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_annotation_file(self):
        '''
        Populate data in "combobox_annotation_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_annotation_file.set('')

        # get the list of the annotation file names
        (_, _, annotation_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the annotation file names in the combobox
        self.combobox_annotation_file['values'] = sorted(annotation_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # clear data in "combobox_annotation_file"
        self.combobox_annotation_file['values'] = []
        self.wrapper_annotation_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        self.populate_combobox_reference_file()

        # load data in "combobox_annotation_file"
        self.populate_combobox_annotation_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_annotation_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_annotation_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTopHatConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_annotation_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the TopHat config file
        if OK:
            message = f'The file {xtophat.get_tophat_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the TopHat config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtophat.create_tophat_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xtophat.create_tophat_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.wrapper_annotation_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the TopHat config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtophat.get_tophat_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtophat.check_tophat_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_tophat_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTopHatConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTransAbyssConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTransAbyssConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_transabyss_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTransAbyssConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTransAbyssConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the Trans-ABySS config file
        if OK:
            message = f'The file {xtransabyss.get_transabyss_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Trans-ABySS config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtransabyss.create_transabyss_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xtransabyss.create_transabyss_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Trans-ABySS config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtransabyss.get_transabyss_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtransabyss.check_transabyss_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_transabyss_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTransAbyssConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTranscriptFilterConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTranscriptFilterConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_transcript_filter_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_rsem_eval_dataset = tkinter.StringVar()
        self.wrapper_rsem_eval_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTranscriptFilterConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_rsem_eval_dataset" and register it with the grid geometry manager
        self.label_rsem_eval_dataset = tkinter.Label(self, text='RSEM-EVAL dataset')
        self.label_rsem_eval_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_rsem_eval_dataset" and register it with the grid geometry manager
        self.combobox_rsem_eval_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_rsem_eval_dataset)
        self.combobox_rsem_eval_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*48)
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_rsem_eval_dataset.bind('<<ComboboxSelected>>', self.combobox_rsem_eval_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_rsem_eval_dataset['values'] = []
        self.wrapper_rsem_eval_dataset.set('')
        self.rsem_eval_dataset_id = None

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_rsem_eval_dataset(self):
        '''
        Populate data in "combobox_rsem_eval_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_rsem_eval_dataset.set('')

        # get the list of the RSEM-EVAL dataset names
        app_list = [xlib.get_rsem_eval_code()]
        (_, _, rsem_eval_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the RSEM-EVAL dataset names in the combobox
        self.combobox_rsem_eval_dataset['values'] = sorted(rsem_eval_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_rsem_eval_dataset"
        self.combobox_rsem_eval_dataset['values'] = []
        self.wrapper_rsem_eval_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_rsem_eval_dataset"
        self.populate_combobox_rsem_eval_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_rsem_eval_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_rsem_eval_dataset" has been selected
        '''

        # get the RSEM-EVAL dataset identification
        (_, _, self.rsem_eval_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_rsem_eval_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTranscriptFilterConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_rsem_eval_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the transcript-filter config file
        if OK:
            message = f'The file {xngshelper.get_transcript_filter_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the transcript-filter config file
        if OK:
            (OK, error_list) = xngshelper.create_transcript_filter_config_file(self.wrapper_experiment_id.get(),  self.rsem_eval_dataset_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the transcript-filter config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xngshelper.get_transcript_filter_config_file())
            self.wait_window(dialog_editor)

            # check the transcript-filter config file
            (OK, error_list) = xngshelper.check_transcript_filter_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_transcript_filter_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTranscriptFilterConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTranscriptomeBlastxConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTranscriptomeBlastxConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_transcriptome_blastx_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_database_dataset = tkinter.StringVar()
        self.wrapper_database_dataset.trace('w', self.check_inputs)
        self.wrapper_protein_database_name = tkinter.StringVar()
        self.wrapper_protein_database_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTranscriptomeBlastxConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_database_dataset" and register it with the grid geometry manager
        self.label_database_dataset = tkinter.Label(self, text='Database dataset')
        self.label_database_dataset.grid(row=1, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_database_dataset" and register it with the grid geometry manager
        self.combobox_database_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_database_dataset)
        self.combobox_database_dataset.grid(row=1, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_protein_database_name" and register it with the grid geometry manager
        self.label_protein_database_name = tkinter.Label(self, text='Protein database')
        self.label_protein_database_name.grid(row=2, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_protein_database_name" and register it with the grid geometry manager
        self.combobox_protein_database_name = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_protein_database_name)
        self.combobox_protein_database_name.grid(row=2, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(40,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(40,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*50)
        self.label_fit.grid(row=6, column=2, padx=(0,0), pady=(40,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=6, column=3, padx=(0,5), pady=(40,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=6, column=4, padx=(5,5), pady=(40,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_database_dataset.bind('<<ComboboxSelected>>', self.combobox_database_dataset_selected_item)
        self.combobox_protein_database_name.bind('<<ComboboxSelected>>', self.combobox_protein_database_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_database_dataset['values'] = []
        self.wrapper_database_dataset.set('')
        self.combobox_protein_database_name['values'] = []
        self.wrapper_protein_database_name.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_database_dataset(self):
        '''
        Populate data in "combobox_database_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_database_dataset.set('')

        # get the list of the database dataset names
        (_, _, database_dataset_name_list) = xdatabase.get_database_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the database dataset names in the combobox
        self.combobox_database_dataset['values'] = sorted(database_dataset_name_list)

    #---------------

    def populate_combobox_protein_database_name(self):
        '''
        Populate data in "combobox_protein_database_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_protein_database_name.set('')

        # get the list of the database file names
        (_, _, database_file_name_list) = xdatabase.get_database_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_database_dataset.get(), file_type='.*phr', passed_connection=True, ssh_client=self.ssh_client)

        # get the list of the database dataset names
        protein_database_name_list = []
        pattern = re.compile('^.*[0-9][0-9]$')
        for database_file_name in database_file_name_list:
            (file_name, _) = os.path.splitext(database_file_name)
            if pattern.match(file_name):
                database_name = file_name[:-3]
            else:
                database_name = file_name
            if database_name not in protein_database_name_list:
                protein_database_name_list.append(database_name)

        # load the database dataset names in the combobox
        self.combobox_protein_database_name['values'] = sorted(protein_database_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_database_dataset"
        self.populate_combobox_database_dataset()

        # clear data in "combobox_protein_database_name"
        self.combobox_protein_database_name['values'] = []
        self.wrapper_protein_database_name.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_database_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_database_dataset" has been selected
        '''

        # load data in "combobox_protein_database_name"
        self.populate_combobox_protein_database_name()

    #---------------

    def combobox_protein_database_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_protein_database_name" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTranscriptomeBlastxConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_database_dataset.get() != '' and self.wrapper_protein_database_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the transcriptome-blastx config file
        if OK:
            message = f'The file {xngshelper.get_transcriptome_blastx_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the transcriptome-blastx config file
        if OK:
            (OK, error_list) = xngshelper.create_transcriptome_blastx_config_file(self.wrapper_database_dataset.get(), self.wrapper_protein_database_name.get(), self.wrapper_experiment_id.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the transcriptome-blastx config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xngshelper.get_transcriptome_blastx_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xngshelper.check_transcriptome_blastx_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_transcriptome_blastx_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTranscriptomeBlastxConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTransrateConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTransrateConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_transrate_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTransrateConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=4, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=4, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=5, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=5, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=5, column=2, columnspan=3, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=6, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=6, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=7, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=7, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=8, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=8, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=9, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=9, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=10, column=0, padx=(15,5), pady=(15,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=10, column=1, padx=(5,5), pady=(15,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=11, column=2, padx=(0,0), pady=(15,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=11, column=3, padx=(0,5), pady=(15,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=11, column=4, padx=(5,5), pady=(15,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')
        self.assembly_dataset_id = None
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('')
        self.combobox_assembly_type['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference dataset names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # clear data in "combobox_assembly_dataset"
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')
        else:
            self.populate_combobox_reference_file()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # load data in "combobox_assembly_dataset"
        self.populate_combobox_assembly_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'readonly'
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTransrateConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != '') and self.wrapper_assembly_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'Invalid pattern.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'A pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the Transrate config file
        if OK:
            message = f'The file {xtransrate.get_transrate_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Transrate config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtransrate.create_transrate_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, selected_file_list, None, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            elif self.read_type == 'PE':
                (OK, error_list) = xtransrate.create_transrate_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list, self.assembly_dataset_id, self.wrapper_assembly_type.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Transrate config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtransrate.get_transrate_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtransrate.check_transrate_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_transrate_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTransrateConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTrimmomaticConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTrimmomaticConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_trimmomatic_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTrimmomaticConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTrimmomaticConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the Trimmomatic config file
        if OK:
            message = f'The file {xtrimmomatic.get_trimmomatic_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Trimmomatic config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtrimmomatic.create_trimmomatic_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xtrimmomatic.create_trimmomatic_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Trimmomatic config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtrimmomatic.get_trimmomatic_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtrimmomatic.check_trimmomatic_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_trimmomatic_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTrimmomaticConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateTrinityConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateTrinityConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_trinity_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_read_dataset = tkinter.StringVar()
        self.wrapper_read_dataset.trace('w', self.check_inputs)
        self.wrapper_file_pattern = tkinter.StringVar()
        self.wrapper_file_pattern.trace('w', self.check_inputs)
        self.wrapper_read_type = tkinter.StringVar()
        self.wrapper_read_type.trace('w', self.check_inputs)
        self.wrapper_specific_chars_1 = tkinter.StringVar()
        self.wrapper_specific_chars_1.trace('w', self.check_inputs)
        self.wrapper_specific_chars_2 = tkinter.StringVar()
        self.wrapper_specific_chars_2.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateTrinityConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_dataset" and register it with the grid geometry manager
        self.label_read_dataset = tkinter.Label(self, text='Read dataset')
        self.label_read_dataset.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_dataset" and register it with the grid geometry manager
        self.combobox_read_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_read_dataset)
        self.combobox_read_dataset.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern" and register it with the grid geometry manager
        self.label_file_pattern = tkinter.Label(self, text='File pattern')
        self.label_file_pattern.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_file_pattern" and register it with the grid geometry manager
        self.entry_file_pattern = tkinter.Entry(self, textvariable=self.wrapper_file_pattern, width=30, validatecommand=self.check_inputs)
        self.entry_file_pattern.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_file_pattern_warning" and register it with the grid geometry manager
        self.label_file_pattern_warning = tkinter.Label(self, text='')
        self.label_file_pattern_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_type" and register it with the grid geometry manager
        self.label_read_type = tkinter.Label(self, text='Read type')
        self.label_read_type.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_type" and register it with the grid geometry manager
        self.combobox_read_type = tkinter.ttk.Combobox(self, width=15, height=4, state='readonly', textvariable=self.wrapper_read_type)
        self.combobox_read_type.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_1" and register it with the grid geometry manager
        self.label_specific_chars_1 = tkinter.Label(self, text='File #1 specific chars')
        self.label_specific_chars_1.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_1" and register it with the grid geometry manager
        self.entry_specific_chars_1 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_1, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_1.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_specific_chars_2" and register it with the grid geometry manager
        self.label_specific_chars_2 = tkinter.Label(self, text='File #2 specific chars')
        self.label_specific_chars_2.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_specific_chars_2" and register it with the grid geometry manager
        self.entry_specific_chars_2 = tkinter.Entry(self, textvariable=self.wrapper_specific_chars_2, width=30, validatecommand=self.check_inputs)
        self.entry_specific_chars_2.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*44)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_read_dataset.bind('<<ComboboxSelected>>', self.combobox_read_dataset_selected_item)
        self.combobox_read_type.bind('<<ComboboxSelected>>', self.combobox_read_type_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')
        self.read_dataset_id = None
        self.wrapper_file_pattern.set('.*fastq')
        self.read_type = None
        self.wrapper_specific_chars_1.set('')
        self.entry_specific_chars_1['state'] = 'disabled'
        self.wrapper_specific_chars_2.set('')
        self.entry_specific_chars_2['state'] = 'disabled'

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_read_type()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identifications list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_read_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_read_dataset(self):
        '''
        Populate data in "combobox_read_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_dataset.set('')

        # get the list of the read dataset names
        (_, _, read_dataset_name_list) = xread.get_read_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the list of the read dataset names in the combobox
        self.combobox_read_dataset['values'] = read_dataset_name_list

    #---------------

    def populate_combobox_read_type(self):
        '''
        Populate data in "combobox_read_type".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_type.set('')

        # load the list of the read dataset names in the combobox
        self.combobox_read_type['values'] =['Single-end', 'Paired-end']

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_read_dataset"
        self.combobox_read_dataset['values'] = []
        self.wrapper_read_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_read_dataset"
        self.populate_combobox_read_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_read_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_dataset" has been selected
        '''

        # get the read dataset identification
        (_, _, self.read_dataset_id) = xread.get_read_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_read_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def combobox_read_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_read_type" has been selected
        '''

        # get the read type code
        if self.wrapper_read_type.get() == 'Single-end':
            self.read_type = 'SE'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.read_type = 'PE'

        # enable or disable the specific chars entries
        if self.wrapper_read_type.get() == 'Single-end':
            self.wrapper_specific_chars_1.set('')
            self.entry_specific_chars_1['state'] = 'disabled'
            self.wrapper_specific_chars_2.set('')
            self.entry_specific_chars_2['state'] = 'disabled'
        elif self.wrapper_read_type.get() == 'Paired-end':
            self.wrapper_specific_chars_1.set('1.fastq')
            self.entry_specific_chars_1['state'] = 'normal'
            self.wrapper_specific_chars_2.set('2.fastq')
            self.entry_specific_chars_2['state'] = 'normal'

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateTrinityConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_file_pattern"
        if not self.check_entry_file_pattern():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_read_dataset.get() != '' and self.wrapper_file_pattern.get() != '' and (self.read_type == 'SE' or self.read_type == 'PE' and self.wrapper_specific_chars_1.get() != '' and  self.wrapper_specific_chars_2.get() != ''):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_file_pattern(self):
        '''
        Check the content of "entry_file_pattern"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_file_pattern" value is a valid pattern of regular expression
        try:
            re.compile(self.wrapper_file_pattern.get())
        except Exception:
            self.label_file_pattern_warning['text'] = 'It is not a valid pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_file_pattern_warning['text'] = 'It is a pattern of regular expression.'
            self.label_file_pattern_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # build the cluster read directory path
        if OK:
            cluster_read_dir = f'{xlib.get_cluster_read_dir()}/{self.wrapper_experiment_id.get()}/{self.read_dataset_id}'

        # get the selected file list
        if OK:
            selected_file_list = []
            command = f'cd {cluster_read_dir}; find . -type f -regex "./{self.wrapper_file_pattern.get()}"'
            (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
            if OK:
                for line in stdout:
                    selected_file_list.append(line.rstrip('\n'))
            else:
                message = f'*** ERROR: Wrong command ---> {command}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            if selected_file_list == []:
                message = f'WARNING: There are not files in the cluster directory {cluster_read_dir} with the pattern {self.wrapper_file_pattern.get()}'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the paired file list when the read type is paired-end
        if OK:
            if self.read_type == 'PE':
                (file_1_list, file_2_list, unpaired_file_list) = xlib.pair_files(selected_file_list, self.wrapper_specific_chars_1.get(), self.wrapper_specific_chars_2.get())
                if unpaired_file_list != []:
                    message = f'ERROR: There are unpaired files: {unpaired_file_list}'
                    tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                    OK = False

        # confirm the creation of the Trinity config file
        if OK:
            message = f'The file {xtrinity.get_trinity_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Trinity config file
        if OK:
            if self.read_type == 'SE':
                (OK, error_list) = xtrinity.create_trinity_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, selected_file_list, None)
            elif self.read_type == 'PE':
                (OK, error_list) = xtrinity.create_trinity_config_file(self.wrapper_experiment_id.get(), self.read_dataset_id, self.read_type, file_1_list, file_2_list)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Trinity config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xtrinity.get_trinity_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xtrinity.check_trinity_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_trinity_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateTrinityConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRecreateVariantCallingConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateVariantCallingConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # assign the text of the "head"
        self.head = f'{xlib.get_variant_calling_name()} - Recreate config file'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_reference_dataset = tkinter.StringVar()
        self.wrapper_reference_dataset.trace('w', self.check_inputs)
        self.wrapper_reference_file = tkinter.StringVar()
        self.wrapper_reference_file.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_assembly_dataset = tkinter.StringVar()
        self.wrapper_assembly_dataset.trace('w', self.check_inputs)
        self.wrapper_assembly_type = tkinter.StringVar()
        self.wrapper_assembly_type.trace('w', self.check_inputs)
        self.wrapper_alignment_dataset = tkinter.StringVar()
        self.wrapper_alignment_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateVariantCallingConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_reference_dataset" and register it with the grid geometry manager
        self.label_reference_dataset = tkinter.Label(self, text='Reference dataset')
        self.label_reference_dataset.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_reference_dataset" and register it with the grid geometry manager
        self.combobox_reference_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_dataset)
        self.combobox_reference_dataset.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_reference_file" and register it with the grid geometry manager
        self.label_reference_file = tkinter.Label(self, text='Reference file')
        self.label_reference_file.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_reference_file" and register it with the grid geometry manager
        self.combobox_reference_file = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_reference_file)
        self.combobox_reference_file.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment id.')
        self.label_experiment_id.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_assembly_dataset" and register it with the grid geometry manager
        self.label_assembly_dataset = tkinter.Label(self, text='Assembly dataset')
        self.label_assembly_dataset.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_assembly_dataset" and register it with the grid geometry manager
        self.combobox_assembly_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_assembly_dataset)
        self.combobox_assembly_dataset.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_assembly_type" and register it with the grid geometry manager
        self.label_assembly_type = tkinter.Label(self, text='Assembly type')
        self.label_assembly_type.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_assembly_type" and register it with the grid geometry manager
        self.combobox_assembly_type = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_assembly_type)
        self.combobox_assembly_type.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_alignment_dataset" and register it with the grid geometry manager
        self.label_alignment_dataset = tkinter.Label(self, text='Alignment dataset')
        self.label_alignment_dataset.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_alignment_dataset" and register it with the grid geometry manager
        self.combobox_alignment_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_alignment_dataset)
        self.combobox_alignment_dataset.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*48)
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(0,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_reference_dataset.bind('<<ComboboxSelected>>', self.combobox_reference_dataset_selected_item)
        self.combobox_reference_file.bind('<<ComboboxSelected>>', self.combobox_reference_file_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_assembly_dataset.bind('<<ComboboxSelected>>', self.combobox_assembly_dataset_selected_item)
        self.combobox_assembly_type.bind('<<ComboboxSelected>>', self.combobox_assembly_type_selected_item)
        self.combobox_alignment_dataset.bind('<<ComboboxSelected>>', self.combobox_alignment_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_reference_dataset['values'] = []
        self.wrapper_reference_dataset.set('')
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_assembly_dataset['values'] = []
        self.wrapper_assembly_dataset.set('NONE')
        self.combobox_assembly_dataset['state'] = 'disabled'
        self.assembly_dataset_id = 'NONE'
        self.combobox_assembly_type['values'] = []
        self.wrapper_assembly_type.set('NONE')
        self.combobox_assembly_type['state'] = 'disabled'
        self.combobox_alignment_dataset['values'] = []
        self.wrapper_alignment_dataset.set('')
        self.alignment_dataset_id = ''

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_reference_dataset(self):
        '''
        Populate data in "combobox_reference_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_dataset.set('')

        # get the list of the reference dataset names
        (_, _, reference_dataset_name_list) = xreference.get_reference_dataset_name_list(self.wrapper_cluster_name.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference dataset names in the combobox
        self.combobox_reference_dataset['values'] = ['NONE'] + sorted(reference_dataset_name_list)

    #---------------

    def populate_combobox_reference_file(self):
        '''
        Populate data in "combobox_reference_file".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_file.set('')

        # get the list of the reference file names
        (_, _, reference_file_name_list) = xreference.get_reference_file_name_list(self.wrapper_cluster_name.get(), self.wrapper_reference_dataset.get(), passed_connection=True, ssh_client=self.ssh_client)

        # load the reference file names in the combobox
        self.combobox_reference_file['values'] = sorted(reference_file_name_list)

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout,_) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_assembly_dataset(self):
        '''
        Populate data in "combobox_assembly_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_assembly_dataset.set('')

        # get the list of the assembly_dataset dataset names
        app_list = [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(), xlib.get_transcript_filter_code(), xlib.get_soapdenovo2_code(), xlib.get_starcode_code()]
        (_, _, assembly_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_assembly_dataset['values'] = sorted(assembly_dataset_name_list)

    #---------------

    def populate_combobox_alignment_dataset(self):
        '''
        Populate data in "combobox_alignment_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_alignment_dataset.set('')

        # get the list of the alignment dataset names
        app_list = [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code(), xlib.get_star_code(), xlib.get_tophat_code()]
        (_, _, alignment_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the alignment dataset names in the combobox
        self.combobox_alignment_dataset['values'] = sorted(alignment_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_reference_dataset"
        self.populate_combobox_reference_dataset()

        # clear data in "combobox_reference_file"
        self.combobox_reference_file['values'] = []
        self.wrapper_reference_file.set('')

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_alignment_dataset"
        self.combobox_alignment_dataset['values'] = []
        self.wrapper_alignment_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_dataset" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_reference_file", "combobox_assembly_dataset" and "combobox_assembly_type"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.combobox_reference_file['values'] = ['NONE']
            self.wrapper_reference_file.set('NONE')

            self.assembly_dataset_id = ''
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = []
            self.wrapper_assembly_dataset.set('')

            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = []
            self.wrapper_assembly_type.set('')
            self.combobox_assembly_type['state'] = 'disabled'

        else:
            self.populate_combobox_reference_file()

            self.assembly_dataset_id = 'NONE'
            self.combobox_assembly_dataset['state'] = 'readonly'
            self.combobox_assembly_dataset['values'] = ['NONE']
            self.wrapper_assembly_dataset.set('NONE')
            self.combobox_assembly_dataset['state'] = 'disabled'

            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_alignment_dataset"
        self.combobox_alignment_dataset['values'] = []
        self.wrapper_alignment_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_reference_file_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_reference_file" has been selected
        '''

        pass

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_assembly_dataset"
        if self.wrapper_reference_dataset.get() == 'NONE':
            self.populate_combobox_assembly_dataset()

        # load data in "combobox_alignment_dataset"
        self.populate_combobox_alignment_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_assembly_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_dataset" has been selected
        '''

        # get the assembly_dataset dataset identification
        (_, _, self.assembly_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_assembly_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

        # load data in "combobox_assembly_type"
        if self.assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or self.assembly_dataset_id.startswith(xlib.get_soapdenovo2_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['CONTIGS', 'SCAFFOLDS']
            self.wrapper_assembly_type.set('')
        elif self.assembly_dataset_id.startswith(xlib.get_transabyss_code()) or self.assembly_dataset_id.startswith(xlib.get_trinity_code()) or self.assembly_dataset_id.startswith(xlib.get_ggtrinity_code()) or self.assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()) or self.assembly_dataset_id.startswith(xlib.get_transcript_filter_code()) or self.assembly_dataset_id.startswith(xlib.get_starcode_code()):
            self.combobox_assembly_type['state'] = 'readonly'
            self.combobox_assembly_type['values'] = ['NONE']
            self.wrapper_assembly_type.set('NONE')
            self.combobox_assembly_type['state'] = 'disabled'

    #---------------

    def combobox_assembly_type_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_assembly_type" has been selected
        '''

        pass

    #---------------

    def combobox_alignment_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_alignment_dataset" has been selected
        '''

        # get the alignment dataset identification
        (_, _, self.alignment_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_alignment_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)


    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateVariantCallingConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_reference_dataset.get() != '' and self.wrapper_reference_file.get() != '' and self.wrapper_experiment_id.get() != '' and self.wrapper_assembly_dataset.get() != '' and self.wrapper_assembly_type.get() != '' and self.wrapper_alignment_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Execute the creation of the config file.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the Variant calling config file
        if OK:
            message = f'The file {xddradseqtools.get_variant_calling_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the Variant calling config file
        if OK:
            (OK, error_list) = xddradseqtools.create_variant_calling_config_file(self.wrapper_experiment_id.get(), self.wrapper_reference_dataset.get(), self.wrapper_reference_file.get(), self.assembly_dataset_id, self.wrapper_assembly_type.get(), self.alignment_dataset_id)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # edit the Variant calling config file
        if OK:

            # edit the config file using "DialogEditor" 
            dialog_editor = gdialogs.DialogEditor(self, xddradseqtools.get_variant_calling_config_file())
            self.wait_window(dialog_editor)

            # check the config file
            (OK, error_list) = xddradseqtools.check_variant_calling_config_file(strict=False)
            if OK:
                message = f'The {xlib.get_variant_calling_name()} config file is OK.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = 'Detected errors:\n\n'
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateVariantCallingConfigFile".
        '''
        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRunBioinfoProcess(tkinter.Frame):

    #---------------

    def __init__(self, main, app):
        '''
        Execute actions correspending to the creation of a "FormRunBioinfoProcess" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container
        self.app = app

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # set the name
        if self.app == xlib.get_bowtie2_code():
            self.name = xlib.get_bowtie2_name()

        elif self.app == xlib.get_busco_code():
            self.name = xlib.get_busco_name()

        elif self.app == xlib.get_cd_hit_est_code():
            self.name = xlib.get_cd_hit_est_name()

        elif self.app == xlib.get_cuffdiff_code():
            self.name = xlib.get_cuffdiff_name()

        elif self.app == xlib.get_cufflinks_cuffmerge_code():
            self.name = xlib.get_cufflinks_cuffmerge_name()

        elif self.app == xlib.get_cuffnorm_code():
            self.name = xlib.get_cuffnorm_name()

        elif self.app == xlib.get_cuffquant_code():
            self.name = xlib.get_cuffquant_name()

        elif self.app == xlib.get_cutadapt_code():
            self.name = xlib.get_cutadapt_name()

        elif self.app == xlib.get_ddradseq_simulation_code():
            self.name = xlib.get_ddradseq_simulation_name()

        elif self.app == xlib.get_express_code():
            self.name = xlib.get_express_name()

        elif self.app == xlib.get_fastqc_code():
            self.name = xlib.get_fastqc_name()

        elif self.app == xlib.get_ggtrinity_code():
            self.name = xlib.get_ggtrinity_name()

        elif self.app == xlib.get_gmap_code():
            self.name = xlib.get_gmap_name()

        elif self.app == xlib.get_gsnap_code():
            self.name = xlib.get_gsnap_name()

        elif self.app == xlib.get_hisat2_code():
            self.name = xlib.get_hisat2_name()

        elif self.app == xlib.get_htseq_count_code():
            self.name = xlib.get_htseq_count_name()

        elif self.app == xlib.get_insilico_read_normalization_code():
            self.name = xlib.get_insilico_read_normalization_name()

        elif self.app == xlib.get_ipyrad_code():
            self.name = xlib.get_ipyrad_name()

        elif self.app == xlib.get_kallisto_code():
            self.name = xlib.get_kallisto_name()

        elif self.app == xlib.get_quast_code():
            self.name = xlib.get_quast_name()

        elif self.app == xlib.get_raddesigner_code():
            self.name = xlib.get_raddesigner_name()

        elif self.app == xlib.get_ref_eval_code():
            self.name = xlib.get_ref_eval_name()

        elif self.app == xlib.get_rnaquast_code():
            self.name = xlib.get_rnaquast_name()

        elif self.app == xlib.get_rsem_eval_code():
            self.name = xlib.get_rsem_eval_name()

        elif self.app == xlib.get_rsitesearch_code():
            self.name = xlib.get_rsitesearch_name()

        elif self.app == xlib.get_soapdenovo2_code():
            self.name = xlib.get_soapdenovo2_name()

        elif self.app == xlib.get_soapdenovotrans_code():
            self.name = xlib.get_soapdenovotrans_name()

        elif self.app == xlib.get_star_code():
            self.name = xlib.get_star_name()

        elif self.app == xlib.get_starcode_code():
            self.name = xlib.get_starcode_name()

        elif self.app == xlib.get_tophat_code():
            self.name = xlib.get_tophat_name()

        elif self.app == xlib.get_transabyss_code():
            self.name = xlib.get_transabyss_name()

        elif self.app == xlib.get_transcript_filter_code():
            self.name = xlib.get_transcript_filter_name()

        elif self.app == xlib.get_transcriptome_blastx_code():
            self.name = xlib.get_transcriptome_blastx_name()

        elif self.app == xlib.get_transrate_code():
            self.name = xlib.get_transrate_name()

        elif self.app == xlib.get_trimmomatic_code():
            self.name = xlib.get_trimmomatic_name()

        elif self.app == xlib.get_trinity_code():
            self.name = xlib.get_trinity_name()

        elif self.app == xlib.get_variant_calling_code():
            self.name = xlib.get_variant_calling_name()

        # assign the text of the "head"
        self.head = f'{self.name} - Run process'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRunBioinfoProcess".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*66)
        self.label_fit.grid(row=1, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=1, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=1, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRunBioinfoProcess" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        if self.app == xlib.get_transcriptome_blastx_code():
            OK = self.check_combobox_cluster_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Run bioinfo process.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the process run
        if OK:
            message = f'The {self.name} process is going to be run.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # execute the process
        if OK:

            # execute the process when it is a Bowtie2 process
            if self.app == xlib.get_bowtie2_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xbowtie2.run_bowtie2_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbowtie2.run_bowtie2_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a BUSCO process
            elif self.app == xlib.get_busco_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xbusco.run_busco_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xbusco.run_busco_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a CD-HIT-EST process
            elif self.app == xlib.get_cd_hit_est_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcdhit.run_cd_hit_est_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcdhit.run_cd_hit_est_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Cuffdiff process
            elif self.app == xlib.get_cuffdiff_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcufflinks.run_cuffdiff_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcufflinks.run_cuffdiff_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Cufflinks-Cuffmerge process
            elif self.app == xlib.get_cufflinks_cuffmerge_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcufflinks.run_cufflinks_cuffmerge_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcufflinks.run_cufflinks_cuffmerge_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Cuffnorm process
            elif self.app == xlib.get_cuffnorm_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcufflinks.run_cuffnorm_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcufflinks.run_cuffnorm_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Cuffquant process
            elif self.app == xlib.get_cuffquant_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcufflinks.run_cuffquant_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcufflinks.run_cuffquant_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a cutadapt process
            elif self.app == xlib.get_cutadapt_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xcutadapt.run_cutadapt_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcutadapt.run_cutadapt_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a ddRADseq simulation process
            elif self.app == xlib.get_ddradseq_simulation_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.run_ddradseq_simulation_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.run_ddradseq_simulation_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a eXpress process
            elif self.app == xlib.get_express_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xexpress.run_express_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xexpress.run_express_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a FastQC process
            elif self.app == xlib.get_fastqc_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xfastqc.run_fastqc_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xfastqc.run_fastqc_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Genome-guided Trinity process
            elif self.app == xlib.get_ggtrinity_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.run_ggtrinity_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.run_ggtrinity_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a GMAP process
            elif self.app == xlib.get_gmap_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xgmap.run_gmap_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xgmap.run_gmap_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a GSNAP process
            elif self.app == xlib.get_gsnap_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xgmap.run_gsnap_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xgmap.run_gsnap_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a HISAT2 process
            elif self.app == xlib.get_hisat2_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xhisat2.run_hisat2_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xhisat2.run_hisat2_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a htseq-count process
            elif self.app == xlib.get_htseq_count_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xhtseq.run_htseq_count_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xhtseq.run_htseq_count_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a insilico_read_normalization process
            elif self.app == xlib.get_insilico_read_normalization_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.run_insilico_read_normalization_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.run_insilico_read_normalization_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a ipyrad process
            elif self.app == xlib.get_ipyrad_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xipyrad.run_ipyrad_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xipyrad.run_ipyrad_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a kallisto process
            elif self.app == xlib.get_kallisto_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xkallisto.run_kallisto_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xkallisto.run_kallisto_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a QUAST process
            elif self.app == xlib.get_quast_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xquast.run_quast_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xquast.run_quast_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a RADdesigner process
            elif self.app == xlib.get_raddesigner_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xraddesigner.run_raddesigner_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xraddesigner.run_raddesigner_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a REF-EVAL process
            elif self.app == xlib.get_ref_eval_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xdetonate.run_ref_eval_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xdetonate.run_ref_eval_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a rnaQUAST process
            elif self.app == xlib.get_rnaquast_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xrnaquast.run_rnaquast_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xrnaquast.run_rnaquast_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a RSEM-EVAL process
            elif self.app == xlib.get_rsem_eval_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xdetonate.run_rsem_eval_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xdetonate.run_rsem_eval_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a rsitesearch process
            elif self.app == xlib.get_rsitesearch_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.run_rsitesearch_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.run_rsitesearch_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a SOAPdenovo2 process
            elif self.app == xlib.get_soapdenovo2_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xsoapdenovo2.run_soapdenovo2_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xsoapdenovo2.run_soapdenovo2_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a SOAPdenovo-Trans process
            elif self.app == xlib.get_soapdenovotrans_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xsoapdenovotrans.run_soapdenovotrans_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xsoapdenovotrans.run_soapdenovotrans_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a STAR process
            elif self.app == xlib.get_star_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xstar.run_star_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xstar.run_star_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a starcode process
            elif self.app == xlib.get_starcode_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xstarcode.run_starcode_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xstarcode.run_starcode_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a TopHat process
            elif self.app == xlib.get_tophat_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtophat.run_tophat_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtophat.run_tophat_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Trans-ABySS process
            elif self.app == xlib.get_transabyss_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtransabyss.run_transabyss_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtransabyss.run_transabyss_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a transcript-filter process
            elif self.app == xlib.get_transcript_filter_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xngshelper.run_transcript_filter_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xngshelper.run_transcript_filter_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a transcriptome_blastx process
            elif self.app == xlib.get_transcriptome_blastx_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xngshelper.run_transcriptome_blastx_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xngshelper.run_transcriptome_blastx_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Transrate process
            elif self.app == xlib.get_transrate_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtransrate.run_transrate_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtransrate.run_transrate_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Trimmomatic process
            elif self.app == xlib.get_trimmomatic_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrimmomatic.run_trimmomatic_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrimmomatic.run_trimmomatic_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Trinity process
            elif self.app == xlib.get_trinity_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.run_trinity_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.run_trinity_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Variant calling process
            elif self.app == xlib.get_variant_calling_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.run_variant_calling_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.run_variant_calling_process, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormRunBioinfoProcess".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRestartBioinfoProcess(tkinter.Frame):

    #---------------

    def __init__(self, main, app):
        '''
        Execute actions correspending to the creation of a "FormRestartBioinfoProcess" instance.
        '''

        # save initial parameters in instance variables
        self.main = main
        self.root = main.root
        self.container = main.container
        self.app = app

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.container)

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # set the name
        if self.app == xlib.get_ddradseq_simulation_code():
            self.name = xlib.get_ddradseq_simulation_name()

        elif self.app == xlib.get_ggtrinity_code():
            self.name = xlib.get_ggtrinity_name()

        elif self.app == xlib.get_insilico_read_normalization_code():
            self.name = xlib.get_insilico_read_normalization_name()

        elif self.app == xlib.get_raddesigner_code():
            self.name = xlib.get_raddesigner_name()

        elif self.app == xlib.get_soapdenovo2_code():
            self.name = xlib.get_soapdenovo2_name()

        elif self.app == xlib.get_soapdenovotrans_code():
            self.name = xlib.get_soapdenovotrans_name()

        elif self.app == xlib.get_trinity_code():
            self.name = xlib.get_trinity_name()

        elif self.app == xlib.get_variant_calling_code():
            self.name = xlib.get_variant_calling_name()

        # assign the text of the "head"
        self.head = f'{self.name} - Restart process'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_experiment_id = tkinter.StringVar()
        self.wrapper_experiment_id.trace('w', self.check_inputs)
        self.wrapper_result_dataset = tkinter.StringVar()
        self.wrapper_result_dataset.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRestartBioinfoProcess".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_experiment_id" and register it with the grid geometry manager
        self.label_experiment_id = tkinter.Label(self, text='Experiment/process')
        self.label_experiment_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_experiment_id" and register it with the grid geometry manager
        self.combobox_experiment_id = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_experiment_id)
        self.combobox_experiment_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_result_dataset" and register it with the grid geometry manager
        self.label_result_dataset = tkinter.Label(self, text='Result dataset')
        self.label_result_dataset.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_result_dataset" and register it with the grid geometry manager
        self.combobox_result_dataset = tkinter.ttk.Combobox(self, width=45, height=4, state='readonly', textvariable=self.wrapper_result_dataset)
        self.combobox_result_dataset.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*46)
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_experiment_id.bind('<<ComboboxSelected>>', self.combobox_experiment_id_selected_item)
        self.combobox_result_dataset.bind('<<ComboboxSelected>>', self.combobox_result_dataset_selected_item)
        self.root.bind('<Return>', self.execute)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_experiment_id['values'] = []
        self.wrapper_experiment_id.set('')
        self.combobox_result_dataset['values'] = []
        self.wrapper_result_dataset.set('')
        self.result_dataset_id = None

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_experiment_id(self):
        '''
        Populate data in "combobox_experiment_id".
        '''

        # clear the value selected in the combobox
        self.wrapper_experiment_id.set('')

        # initialize the experiment identification list
        experiment_id_list = []

        # get the experiment identifications
        command = f'ls {xlib.get_cluster_result_dir()}'
        (OK, stdout, _) = xssh.execute_cluster_command(self.ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    experiment_id_list.append(line)

        # check if there are any experimment identifications
        if experiment_id_list == []:
            message = f'The cluster {self.wrapper_cluster_name.get()} does not have experiment data.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the experiment identifications in the combobox
        if self.app in [xlib.get_ddradseq_simulation_code()]:
            if xlib.get_design_dataset_name() in experiment_id_list:
                experiment_id_list = [xlib.get_design_dataset_name()]
            else:
                experiment_id_list = []
        else:
            if xlib.get_design_dataset_name() in experiment_id_list:
                experiment_id_list.remove(xlib.get_design_dataset_name())
        self.combobox_experiment_id['values'] = sorted(experiment_id_list)

    #---------------

    def populate_combobox_result_dataset(self):
        '''
        Populate data in "combobox_result_dataset".
        '''

        # clear the value selected in the combobox
        self.wrapper_result_dataset.set('')

        # get the list of the assembly dataset names
        app_list = [self.app]
        (_, _, result_dataset_name_list) = xresult.get_result_dataset_name_list(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), 'uncompressed', app_list, passed_connection=True, ssh_client=self.ssh_client)

        # load the assembly dataset names in the combobox
        self.combobox_result_dataset['values'] = sorted(result_dataset_name_list)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_experiment_id"
        self.populate_combobox_experiment_id()

        # clear data in "combobox_result_dataset"
        self.combobox_result_dataset['values'] = []
        self.wrapper_result_dataset.set('')

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_experiment_id_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_experiment_id" has been selected
        '''

        # set cursor to show busy status
        self.root.config(cursor='watch')
        self.root.update()

        # load data in "combobox_result_dataset"
        self.populate_combobox_result_dataset()

        # set cursor to show normal status
        self.root.config(cursor='')
        self.root.update()

    #---------------

    def combobox_result_dataset_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_result_dataset" has been selected
        '''

        # get the assembly dataset identification
        (_, _, self.result_dataset_id) = xresult.get_result_dataset_id(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.wrapper_result_dataset.get(), status='uncompressed', passed_connection=True, ssh_client=self.ssh_client)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRestartBioinfoProcess" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_experiment_id.get() != ''  and self.wrapper_result_dataset.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self, event=None):
        '''
        Restart a bioinformatic process.
        '''

        # if "button_execute" is disabled, exit function
        if str(self.button_execute['state']) == 'disabled':
            return

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the process run
        if OK:
            message = f'The {self.name} process is going to be run.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # execute the process
        if OK:

            # execute the process when it is a ddRADseq simulation process
            if self.app == xlib.get_ddradseq_simulation_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.restart_ddradseq_simulation_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.restart_ddradseq_simulation_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Genome-guided Trinity process
            elif self.app == xlib.get_ggtrinity_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.restart_ggtrinity_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.restart_ggtrinity_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a insilico_read_normalization process
            elif self.app == xlib.get_insilico_read_normalization_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.restart_insilico_read_normalization_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.restart_insilico_read_normalization_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a RADdesigner process
            elif self.app == xlib.get_raddesigner_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xraddesigner.restart_raddesigner_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xraddesigner.restart_raddesigner_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a SOAPdenovo2s process
            elif self.app == xlib.get_soapdenovo2_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xsoapdenovo2.restart_soapdenovo2_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xsoapdenovo2.restart_soapdenovo2_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a SOAPdenovo-Trans process
            elif self.app == xlib.get_soapdenovotrans_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xsoapdenovotrans.restart_soapdenovotrans_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xsoapdenovotrans.restart_soapdenovotrans_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Trinity process
            elif self.app == xlib.get_trinity_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xtrinity.restart_trinity_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xtrinity.restart_trinity_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

            # execute the process when it is a Variant calling process
            elif self.app == xlib.get_variant_calling_code():
                dialog_log = gdialogs.DialogLog(self, self.head, xddradseqtools.restart_variant_calling_process.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xddradseqtools.restart_variant_calling_process, args=(self.wrapper_cluster_name.get(), self.wrapper_experiment_id.get(), self.result_dataset_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormRestartBioinfoProcess".
        '''

        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the classes related to BioInfo application forms in gui mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
