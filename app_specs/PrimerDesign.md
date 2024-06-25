
# Application specification: PrimerDesign

This is the application specification for service with identifier PrimerDesign.

The backend script implementing the application is [App-PrimerDesign.pl](../service-scripts/App-PrimerDesign.pl).

The raw JSON file for this specification is [PrimerDesign.json](PrimerDesign.json).

This service performs the following task:   Use Primer3 to design primers to given sequence

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| input_type | How to interpret sequence_input | ARRAY(0x56109dc1b5d0)  | :heavy_check_mark: |  |
| sequence_input | DNA sequence data | string  | :heavy_check_mark: |  |
| SEQUENCE_ID | Sequence ID | string  |  |  |
| SEQUENCE_TARGET | Amplification target | array  |  |  |
| SEQUENCE_INCLUDED_REGION | Region where primers can be picked | array  |  |  |
| SEQUENCE_EXCLUDED_REGION | Region where primers cannot overlap | array  |  |  |
| SEQUENCE_OVERLAP_JUNCTION_LIST | Start position and length of region that primers must flank | array  |  |  |
| PRIMER_PICK_INTERNAL_OLIGO | Pick an internal oligo (hybridization probe). | number  |  |  |
| PRIMER_PRODUCT_SIZE_RANGE | Min, max product size | array  |  |  |
| PRIMER_NUM_RETURN | Max num primer pairs to report | integer  |  |  |
| PRIMER_MIN_SIZE | Min primer length | integer  |  |  |
| PRIMER_OPT_SIZE | Optimal primer length | integer  |  |  |
| PRIMER_MAX_SIZE | Maximum primer length | integer  |  |  |
| PRIMER_MAX_TM | Maximum primer melting temperature | number  |  |  |
| PRIMER_MIN_TM | Minimum primer melting temperature | number  |  |  |
| PRIMER_OPT_TM | Optimal primer melting temperature | number  |  |  |
| PRIMER_PAIR_MAX_DIFF_TM | Max Tm difference of paired primers | number  |  |  |
| PRIMER_MAX_GC | Maximum primer GC percentage | number  |  |  |
| PRIMER_MIN_GC | Minimum primer GC percentage | number  |  |  |
| PRIMER_OPT_GC | Optimal primer GC percentage | number  |  |  |
| PRIMER_SALT_MONOVALENT | Concentration of monovalent cations (mM) | number  |  |  |
| PRIMER_SALT_DIVALENT | Concentration of divalent cations (mM) | number  |  |  |
| PRIMER_DNA_CONC | Annealing Oligo Concentration (nM) | number  |  |  |
| PRIMER_DNTP_CONC | Concentration of dNTPs (mM) | number  |  |  |

