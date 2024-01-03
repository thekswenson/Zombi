

<img src="https://github.com/AADavin/Zombi/blob/master/Images/ZombiLogo.png" alt="zombilogo" height = "150" width="300"/>

### ** Zombi-fork + ZombiPhy: A more flexible version of Zombi with Simphy integration **

----------

### **Introduction** ###

**Fork Description** This fork is an effort by the Mirarab lab and Dr. Krister M. Swenson to create a simulation program which is able to simulate the evolution of large-scale eukaryotic genomes from a pre-set ancestral root genome.

**Zombi** Zombi is a multilevel simulation program that allows for the simulation of species trees, gene/intergene trees, and sequence evolution. Importantly, the gene/intergene tree simulation incorporates large scale genomic rearrangements into its simulation, such as transpositions and inversions. 

**SimPhy** SimPhy is a fast, open source simulation program that can simulate multiple gene families evolving under gene duplication and loss, gene conversion, and most importantly, incomplete lineage sorting. 

**ZombiPhy** 

----------

### **Installation** ###

**Installing Zombi**
First, clone the repository to your computer using git clone. If you are unfamiliar with git, please see the 
github docs here https://docs.github.com/en/get-started/using-git

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

Use mamba to install the following packages (Note too self, replace with .yaml file):
    
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
        
**Installing Simphy**
There are two methods for installing Simphy, Compilation and pre-compled binaries. We reccomend using the pre-compiled binaries, as they are by far the simplest method for installing SimPhy. You can find the precompiled binaries here: https://github.com/adamallo/SimPhy/releases/latest 

If for whatever reason you wish to compile SimPhy, please follow the instructions in the SimPhy manual here: https://github.com/adamallo/SimPhy/wiki/Manual#4-obtaining-simphy

## **Running ZombiPhy** ##

This is a guide for how to run the the ZombiPhy 

### **Running Zombi the normal way** ###

**There is a detailed Wiki that explains how Zombi works [here](https://github.com/AADavin/ZOMBI/wiki)** and it takes around 15 minutes of your time reading it! But if you want to launch it
right away, just read this.  

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

### Limitations of ZombiPhy future tasks ###

- The output of ZombiPhy is not a valid input for Zombi's sequence simulation mode (S). This will need to be a major task for future collaborators. 
- The pseudogeneization option for Zombi, and by extension ZombiPhy, is broken. For now, this just means that the pseudogeneization option in the ZombiPhy/Zombi genome parameters file must always be left at 0.
- Inputting locus trees with horizontal transfers into SimPhy is broken. For now, this means that the transfer rate in the ZombiPhy/Zombi genome parameters file must always be left at 0.
- Inputting locus trees with losses into SimPhy is broken. This only causes problems in a rare circumstance where a duplication occurs but one of the duplicates is lost. Under these circumstances, the branch length of the surviving duplicate will be inaccurate. For now, this inaccuracy is ignored, but it must be addressed in the future. 
- Zombi cannot support full eukaryotic genomes. The internal infrastructure needed for inputting in an entire, multi-chromosome genome into Zombi simply does not exist. This will have to be implemented in a future version of the program. 
    
    
----------

Writen by Keegan R. Flanagan

Please, if you have any doubts or need a hand, **just contact me here: keflanagan@ucsd.edu**

----------

### **References** ###

**Zombi: A phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.**

Adrian A. Davin, Theo Tricou, Eric Tannier, Damien M. de Vienne, Gergely J. Szollosi

Bioinformatics, btz710, https://doi.org/10.1093/bioinformatics/btz710

**SimPhy: phylogenomic simulation of gene, locus, and species trees. Systematic biology.**

Diego Mallo, Leonardo de Oliveira Martins, and David Posada

Systematic biology

----------



