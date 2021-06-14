Toward a methodology for evaluating DNA variants in nuclear families
================================================================================================================
The [code_used_to_analyze_data.pdf](https://github.com/dmiller903/PedFam/blob/master/code_used_to_analyze_data.pdf)
file contained in this repository provides the commands and code used to
process the VCF files analyzed as part of the manuscript, 
"Toward a methodology for evaluating DNA variants in nuclear families"

The scripts contained within this repository were executed with the docker container
[compound-het-vip](https://hub.docker.com/r/dmill903/compound-het-vip).
All scripts except for `create_gnomAD_1K_cadd_file.py` and  `keep_passed_variants.py` were adapted from
[CompoundHetVIP](https://github.com/dmiller903/CompoundHetVIP). The
example
[here](https://github.com/dmiller903/CompoundHetVIP/blob/master/CompoundHetVIP_example.pdf)
explains, in detail, how the adapted scripts are used. The publication for
CompoundHetVIP can be found here:
[https://doi.org/10.12688/f1000research.26848.2](https://f1000research.com/articles/9-1211)