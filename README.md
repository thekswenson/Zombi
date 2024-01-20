
![ZombiPhy_logo](https://github.com/thekswenson/Zombi/blob/Simphy_integration/Images/ZombiPhy_logo.jpg)

### ** Zombi-fork + ZombiPhy: A more flexible version of Zombi with SimPhy integration **

----------

### **Introduction** ###

**Fork Description:**

This fork is an effort by the Mirarab lab and Dr. Krister M. Swenson to create a simulation program which is able to simulate the evolution of large-scale eukaryotic genomes from a pre-set ancestral root genome.

**Zombi:** 

Zombi is a multilevel simulation program that allows for the simulation of species trees, gene/intergene trees, and sequence evolution. Importantly, the gene/intergene tree simulation incorporates large scale genomic rearrangements into its simulation, such as transpositions and inversions. 

**SimPhy:** 

SimPhy is a fast, open source simulation program that can simulate multiple gene families evolving under gene duplication and loss, gene conversion, and incomplete lineage sorting. 

**ZombiPhy:** 

ZombiPhy is a Python based pipeline that integrates the outputs of Zombi and SimPhy together in order to create a more complete simulation that considers both large scale genomic rearrangements and incomplete lineage sorting. 

----------

### **Installation** ###

**Installing Zombi:**

First, clone the repository to your computer using git clone. If you are unfamiliar with git, please see the 
github docs here https://docs.github.com/en/get-started/using-git.

    git clone https://github.com/thekswenson/Zombi.git

Second, we are going to create a conda environment in which to run Zombi. Here we are using a C++ Reimplementation 
of conda called mamba due to its increased speed and stability. You can learn more about mamba and how to 
install it here: https://mamba.readthedocs.io/en/latest/index.html. If you are more comfortable using 
miniconda or anaconda, please feel free to replace any instance of "mamba" in the commands below with "conda", 
and it will work just fine.

    mamba create -n zombi
    mamba activate zombi

Then we are going to add a couple of additional channels to our new environment. We are using the conda
command for this because mamba cannot alter the configuration. 

    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels nodefaults

Use mamba to install the following packages:
    
    mamba install python ete3 networkx flask gffutils

The packages BCBio and Pyvolve are also both required for zombi, but bioconda does not currently support
the most up to date version of these packages. Therefore, we will download them with pip. 

    pip install bcbio-gff
    pip install pyvolve

Finally, in order to install the Zombi package itself, cd into the Zombi directory downloaded from git
and use pip to install Zombi. The contents of the setup.cfg and pyproject.toml files will automatically be 
used to configure your environment.

    cd Zombi
    pip install -e .
        
**Installing SimPhy:**

There are two methods for installing SimPhy, Compilation and pre-compiled binaries. We reccommend using the pre-compiled binaries, as they are by far the simplest method for installing SimPhy. You can find the precompiled binaries here: https://github.com/adamallo/SimPhy/releases/latest 

If you wish to compile SimPhy, please follow the instructions in the SimPhy manual here: https://github.com/adamallo/SimPhy/wiki/Manual#4-obtaining-simphy

## **Running ZombiPhy** ##

**Overview:**

The ZombiPhy pipeline can be broken down into 5 main steps:
1. SimPhy generates a species tree.
2. The species tree is converted into a format identical to the output of Zombi's **T** mode.
3. Zombi is run in **G** mode using the species tree as input. 
4. Gene trees generated from step 4 are fed into SimPhy as locus trees. This generates new gene trees that incorporate incomplete lineage sorting.
5. The new, SimPhy generated gene trees are fed back into Zombi which is then run in its **S** mode, which generates sequences from the tree. 

The first 4 steps are done using the **ZombiPhy.py** script, which automates the steps for ease of use. The final step is done using a seperate Python script, **Name Pending**, as it is **very computationally expensive**

**Parameters:**

A majority of the parameters fed to ZombiPhy are specified in parameter tsv files rather than in arguments. ZombiPhy parameter files, such as the one in the Parameters folder, are split into two section. 

The first section is a set of arguments to SimPhy. This allows you to change parameters such as the number of replicates and the distribution from which to draw the gene/species specific rate heterogeneity parameters. Note that you should be able to add any of the arguments specified in https://github.com/adamallo/SimPhy/wiki/Manual to this section of the parameters file to add them to the SimPhy call, though this has not been thouroughly tested. For more information on using SimPhy, please see https://github.com/adamallo/SimPhy/wiki/Manual. 

The second section details the parameters that will be inputted into Zombi during genome simulation. This allows you to change parameters such as the duplication, loss, inversion, and transposition rates. Note that unlike the SimPhy parameters section, you cannot add nor remove any of the Zombi parameters. For more information on these parameters, please see the genome generation section of the Zombi wiki at https://github.com/AADavin/Zombi/wiki/Generating-Genomes. 

**Running ZombiPhy:**

ZombiPhy accepts 5 different arguments:
1. -s_loc is the location of the precompiled binaries for SimPhy (if you compiled SimPhy, then you can replace this with the command used to call SimPhy). 
2. -params is the location of the zombiphy parameters file. 
3. -o is the location in which to save the output files. 
4. -g_loc is the location of the root genome gff file. In this example, 
5. -fc is the feature choice. You can choose to look at either individual CDS regions or entire genes.

Here is an example run of ZombiPhy using the default parameters in Parameters/ and the genome data from chicken chromosome 4 (Available here: https://ftp.ensembl.org/pub/release-110/gff3/gallus_gallus/)

    Python ZombiPhy.py -s_loc ../SimPhy_1.0.2/bin/simphy_mac64 -params Parameters/ZombiPhyParameters.tsv -o ../output -g_loc ../genomes/chicken/chr4.gff -fc gene

**Running name pending** ... 

### Running Zombi Independently ###

**There is a detailed Wiki that explains how Zombi works [here](https://github.com/AADavin/ZOMBI/wiki)** and it takes around 15 minutes of your time reading it! But if you want to launch it right away, just read this.  

There are **three main modes** to run Zombi: **T** (species Tree), **G** (Genomes) and  **S** (Sequences) 

You must run the computations in sequential order. This means that:
Computing **genomes** requires having computed previously a species tree computed with the T mode. 
Computing **sequences** requires having computed previously **genomes** with the G mode.

<img src="https://github.com/AADavin/Zombi/blob/master/Images/ZombiGraphicAbstract.png" alt="zombipipeline" height = "300" width="800"/>

The parameters are read from a .tsv file that can be modified with any text editor. 

To start using Zombi just write:

    Zombi T ./Parameters/SpeciesTreeParameters.tsv ../Output_folder

These will generate a Species Tree along with some other files in ./Output folder/T

Then, you can simulate the evolution of genomes in that species tree using:

    Zombi G ./Parameters/GenomeParameters.tsv ../Output_folder

Make sure that the Output_folder is the same that you were using when you generated the Species Tree! 
Zombi will simulate the evolution of genomes inside the whole species tree and it will print a very detailed
information about them in ./Output folder/G

Then, you can simulate the evolution of sequences for each gene in that species tree using:

    Zombi S ./Parameters/SequenceParameters.tsv ../Output_folder

### Example SimPhy Run ###

There are many different ways in which to run SimPhy. We have chosen to include two examples here which illustrate how it is used within ZombiPhy. For a more detailed description of how to run SimPhy, please see https://github.com/adamallo/SimPhy/wiki/Manual. 

Below is what can be seen as a traditional run of SimPhy. This will generate a species tree, a set of 1000 locus trees, and a set of 1000 gene trees, for 3 replicates: 

    SimPhy_1.0.2/bin/simphy_mac64 -SL f:10 -SP f:20000 -SB f:0.0001 -SU f:0.001 -SG f:1 -HS LN:1.5,1 -HL LN:1.5,1 -HG LN:1.5,1 -RS 3 -RL 1000 -CS 220 -V 3 -O output/simphy_output/species_tree

1. -SL is the number of leaves for the species tree.
2. -SP is the tree wide population size.
3. -SB is the speciation rate.
4. -SU is the substitution rate.
5. -SG is the tree wide generation time.
6. -HS, -HL, -HG are the species, lineage, and gene specific rate heterogeneity modifiers respectively.
7. -RS is the number of replicates.
8. -RL is the number of locus trees.
9. -CS is the random seed.
10. -V is the verbosity.
11. -O is the output location.

When this code is run within ZombiPhy, -RL is set to 1 as we will be using Zombi to generate the initial gene trees. 

ZombiPhy runs SimPhy twice. The second run uses the -LR argument to input a set of locus trees created by Zombi into simphy as a nexus file. An example of this has not been included as it is difficult to run SimPhy in -LR mode without running the whole pipeline. 

### Limitations of ZombiPhy and future tasks ###

- The pseudogeneization option for Zombi, and by extension ZombiPhy, is broken. For now, this just means that the pseudogeneization option in the ZombiPhy/Zombi genome parameters file must always be left at 0.
- Inputting locus trees with horizontal transfers into SimPhy is broken. For now, this means that the transfer rate in the ZombiPhy/Zombi genome parameters file must always be left at 0.
- Zombi cannot support full eukaryotic genomes. The internal infrastructure needed for inputting in an entire, multi-chromosome genome into Zombi does not yet exist. This will have to be implemented in a future version of the program.
- There appears to be a bug in Zombi where certain events affect far too many genes (e.g. an inversion happening to 95% of the genes). This is likely due to a problem with how Zombi chooses the number of genes to be affected by each individual event.  
    
    
----------

Writen by Keegan R. Flanagan

Please, if you have any doubts or need a hand, **please contact me here: keflanagan@ucsd.edu**

----------

### **References** ###

**Zombi: A phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.**

Adrian A. Davin, Theo Tricou, Eric Tannier, Damien M. de Vienne, Gergely J. Szollosi

Bioinformatics, btz710, https://doi.org/10.1093/bioinformatics/btz710

**SimPhy: phylogenomic simulation of gene, locus, and species trees. Systematic biology.**

Diego Mallo, Leonardo de Oliveira Martins, and David Posada

Systematic biology

----------



