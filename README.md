# NGScloud2

NGScloud2, a major version of NGScloud, is a bioinformatic system developed to analyze data of
NGS experiments using the cloud computing services of Amazon - Elastic Compute Cloud (EC2)-
that permit the access to ad hoc computing infrastructure scaled according to the complexity of
each experiment, so its costs and times can be optimized. The application provides a user-friendly
front-end to easily operate Amazon's hardware resources and to control workflows of bioinformatic
analysis of *de novo* RNA-seq, reference-based RNA-seq and functional annotation specifically oriented
to plant species. This last workflow encapsulates our standalone application TOA (https://github.com/GGFHF/TOA/)
so it can run in EC2. NScloud2 also has technical improvements as the usage of spot instance
(they can lead to significant cost savings) or the last instances types (M5, C5 and R5).

Refer to the NGScloud2-manual.pdf in the "Package" folder for installation instructions, a description of
the software and examples of use (https://github.com/GGFHF/NGScloud2/blob/master/Package/NGScloud2-manual.pdf).
Also, you can see the paper:

Mora-Márquez F, Vázquez-Poletti JL, López de Heredia U (2018).
RNA-seq analysis of non-model species using cloud computing. *Bioinformatics, 34*(19), 3405–3407.
DOI: https://doi.org/10.1093/bioinformatics/bty363

This software has been developed by:

    GI Genética, Fisiología e Historia Forestal
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politécnica de Madrid
    
    https://github.com/ggfhf/

### Disclaimer

The software package NGScloud2 is available for free download from the GitHub software repository
(https://github.com/GGFHF/NGScloud2) under GNU General Public License v3.0.

The usage of the services of Amazon Web Services, such as the Elastic Compute Cloud (EC2) and
the Elastic Block Store (EBS) that uses NGScloud2, entails expenses. You can inform yourself on
the Amazon website: htpps://aws.amazon.com. The NGScloud2 development team does not assume any
responsibility or liability for the expenses derived from the utilization of the services provided
by Amazon Web Services by the NGScloud2 users.
