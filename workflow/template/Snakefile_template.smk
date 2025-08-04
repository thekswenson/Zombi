# coding: utf-8
"""
This Snakefile is used to run Zombi simulations over a range of parameters.

It is an empty workflow that includes the rules installed in the mamba
environment. The configfile is loaded and then the `all` rule fills in the
target files using `buildAllTargetList` from ZOMBI_SNAKEFILE.
"""
from zombi.snakemake.parameters import ZOMBI_SNAKEFILE

configfile: 'config.yaml'
include: ZOMBI_SNAKEFILE


#-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-
# Targets
#-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-

rule all:
  input:
    buildAllTargetList