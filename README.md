# PCR Primer Design Script

This script designs PCR primers for specified genomic regions using **Primer3** and a **local human reference genome FASTA**.

## Overview
- Accepts a CSV file of genomic target regions
- Extracts flanking sequence from a local FASTA using `pyfaidx`
- Designs primers with **strict Primer3 parameters**, if no primers can be designed according to the given parameters, a second set of relaxed parameters are tried.
- Outputs primers **with and without P5/P7 adapter flaps**
- Reports primer metrics and the full amplicon sequence

## Requirements
- Python 3
- `primer3-py==2.2.0`
- `pyfaidx==0.9.0.3`
- `pandas`
- Local genome FASTA (e.g. `hg38.fa`)

## Input
CSV file with the following columns:

* site_name : an identifier for the site to be amplified
* chrom : chromosome name
* start : start coordinate of target withing chromosome of interest
* end : end coordinate of target within chromosome of interest

The "target" region is defined by start and end coordinates - this is the minimal sequence you want to amplify and primers will be designed outside of the target region.

## Parameters
Parameters can be set within the header of the python script 'design_primer_panel.py'

* `GENOME_FASTA`: path to human genome fasta file
* `FLANK_SIZE`: number of bases on either side of target region to retrieve from human genome sequence - primers will be designed somewhere within these flanking regions (default=200)
* `PRODUCT_SIZE_RANGE`: minimum and maximum possible lengths of the PCR amplicon to be designed (without P5 and P7 flaps included) (defulault [150, 250])
* `PRIMER_3_PARAMS_STRICT`: ideal parameters for Primer3 primer designs. These parameters will be tried first. Default values are shown below. If no Primer3 fails to design primers with these parameters,     relaxed parameters will be used.
  
   ```
   PRIMER3_PARAMS_STRICT = {
     "PRIMER_MIN_SIZE": 18,
     "PRIMER_OPT_SIZE": 23,
     "PRIMER_MAX_SIZE": 26,
  
     "PRIMER_MIN_TM": 63.0,
     "PRIMER_OPT_TM": 65.0,
     "PRIMER_MAX_TM": 70.0,
  
     "PRIMER_MIN_GC": 40.0,
     "PRIMER_MAX_GC": 60.0,
  
     "PRIMER_NUM_RETURN": 1,
  }
  ```

* `PRIMER_3_PARAMS_RELAXED`: Fallback parameters to be used by Primer3 if fails to design primers with optimal parameters. Default values shown below.

  ```
  PRIMER3_PARAMS_RELAXED = {
    **PRIMER3_PARAMS_STRICT,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 65.0,
    "PRIMER_MAX_TM": 70.0,
  
    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,
  
    "PRIMER_MAX_SIZE": 30
  }
  ```
* `P5_FLAP` : P5 flap to append to forward primer 
* `P7_FLAP` : P7 flap tp append to reverse primer 

## Output
CSV file containing, for each site:
- Primer sequences (with and without P5/P7 flaps)
- Primer melting temperatures and GC content
- PCR product size
- Amplicon genomic coordinates and length
- Amplicon sequence
- Design status:
  - `OK` – strict parameters succeeded
  - `RELAXED` – relaxed parameters used
  - `NO_PRIMERS` – no primers found
  - `ERROR` – runtime error for that site

## Usage
```python design_primer_panel.py input.csv output.csv```

Examples of input files and output files are included for a set of off target sites for a CYBB-targeting gRNA ('CYBB_input.csv', 'CYBB_output.csv')
