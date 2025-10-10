"""
Export Zombi simulations to MCGP input formats.
"""

rule export_Zombi_to_FFGC:
  """
  Create FFGC input files in `ffgc_input/` in the simulation directory.
  """
  input:
    SIMDIR + '/sequences/{zparams}/S/Genes',
    SIMDIR + '/sequences/{zparams}/G',  #/Genomes',

  output:
    directory(SIMDIR + '/sequences/{zparams}/ffgc_input')

  shell:
    'zombiExporter ffgc ' + SIMDIR + '/sequences/{wildcards.zparams} {output}'
