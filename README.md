# Primer Design Service

## Overview

The Primer Design Service utilizes Primer3[1-5] to design primers from a given input sequence under a variety of temperature, size, and concentration constraints. 
Primer3 can be used to design primers for several uses including, PCR (Polymerase Chain Reaction) primers, hybridization probes, and sequencing primers. 
The service accepts a nucleotide sequence (pasted in, or select a FASTA formatted file from the workspace) and allows users to specify design. 
Several options enable deliniating which regions of the sequence to search for optimal primers such as: permited ranges, forbinden ranges, and a target region which primers are required to flank.
After specifying an appropriate output folder and clicking “submit”, the primer design is queued as a “job” to process in the Jobs information box on the bottom right of the page. 
Once the job has successfully completed, the output file will appear in the workspace, allowing the user to choose from a list of appropriate primer pairs qualified by scores on several optimality criteria. 



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [PrimerDesign](app_specs/PrimerDesign.md)


## See also

* [Primer Design Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/primer_design_service.html)
* [Primer Design Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/PrimerDesign.html)
* [Primer Design Service Tutorial](https://www.bv-brc.org/docs//tutorial/primer_design/primer_design.html)



## References


1.	Rozen S, Skaletsky H (2000) Primer3 on the WWW for general users and for biologist programmers. Methods Mol Biol 132(3):365–386
2.	Untergasser A, Nijveen H, Rao X, Bisseling T, Geurts R, Leunissen JA (2007) Primer3Plus, an enhanced web interface to Primer3. Nucleic Acids Res 35(Web Server issue):71–74
3.	You FM, Huo N, Gu YQ, Luo MC, Ma Y, Hane D, Lazo GR, Dvorak J, Anderson OD (2008) BatchPrimer3: a high throughput web application for PCR and sequencing primer design. BMC Bioinformatics 9:253
4.	https://github.com/primer3-org
5.	https://primer3.org/manual.html#PRIMER_DNA_CONC
