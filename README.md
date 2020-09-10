# OligoMiner

[![DOI](https://zenodo.org/badge/DOI/10.1073/pnas.1714530115.svg)](http://dx.doi.org/10.1073/pnas.1714530115)

## Overview

This repository contains the code for the [OligoMiner](http://dx.doi.org/10.1073/pnas.1714530115) tool.

If you are looking to use probe sequences that we have already generated for various genome assemblies (hg19, hg38, mm9, mm10, dm3, dm6, ce6, ce11, danRer10, tair10), you can download those on our [website](http://genetics.med.harvard.edu/oligopaints). If you would like to run the OligoMiner tool yourself, please see below for instructions.


## Installing OligoMiner dependencies

1. Make sure you have [conda](https://docs.conda.io/en/latest/miniconda.html) installed. 

2. Clone this repo, then create and activate the provided [environment](./environment.yml):

```
$ git clone https://github.com/beliveau-lab/OligoMiner.git
$ cd OligoMiner
$ conda env create -f environment.yml
$ conda activate probeMining
```

This will install the following packages and their dependencies:

* [Python 2.7.15](https://www.python.org/downloads/release/python-2715/)
* [Biopython](https://biopython.org/)
* [scikit-learn](https://scikit-learn.org/stable/)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [JELLYFISH](https://www.cbcb.umd.edu/software/jellyfish/)
* [NUPACK](http://www.nupack.org/)

### Note about operating systems

OligoMiner is a set of command-line scripts developed on Python 2.7 that can easily be executed from a [Bash Shell](https://en.wikipedia.org/wiki/Bash_(Unix_shell)). If you are using standard Linux or Mac OS X sytsems, we expect these instructions to work for you.

If you are using Windows 10, we recommend enabling [Ubuntu on Windows 10](https://ubuntu.com/tutorials/ubuntu-on-windows), a full Linux distribution, and then running OligoMiner in the Ubuntu terminal.

## Running OligoMiner locally

To make sure all of your dependencies are set up properly, below we will run you through the pipeline using some small example datasets.

### Running scripts on the example files

1. To run the `blockParse.py` script on a .fa file, you can run the following command:

		python blockParse.py -f 3.fa

	This produces a .fastq file (`3.fastq`) containing all identified probe sequences matching your provided criteria. To see additional command line arguments available for this script, you can run the python file with the `-h` argument (i.e. `python blockParse.py -h').

2. NGS alignment. For example, you can use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the newly generated set of candidate probes by running:

		bowtie2 -x /path_to_hg38_index/hg38 -U 3.fastq --no-hd -t -k 100 --very-sensitive-local -S 3_u.sam

	or

		bowtie2 -x /path_to_hg38_index/hg38 -U 3.fastq --no-hd -t -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S 3.sam

	... where "path_to_hg38_index" is replaced with the path to the bowtie2 indices for your genome of interest. These commands produce .sam files (`3_u.sam` and `3.sam`) containing sequence alignment information, but require genome builds as described in the previous section. If you are just testing your scripts to make sure they are working properly, we have already provided the output `3_u.sam` and `3.sam` files in the example files directory for you to use to test subsequent scripts.

3. To process the .sam file produced by sequence alignment, use the `outputClean.py` script:

		python outputClean.py -u -f 3_u.sam

	or, optionally (requires sklearn for the LDA model, see above)

		python outputClean.py -T 42 -f 3.sam

	13 of 13 of the candidate probes should pass the first command (and 12 of 13 candidate probes should pass the specificity filtering with the 42C LDA model in the second command). To see additional command line arguments available for this script, you can run the python file with the `-h` argument (i.e. `python outputClean.py -h').

4. [Optional] Now, you can use `kmerFilter.py` to screen your probes against high abundance kmers (requires [Jellyfish](http://www.genome.umd.edu/jellyfish.html) to be installed and in your path, and a Jellyfish dictionary, see instructions above).

		python kmerFilter.py -f 3_probes.bed -m 18 -j sp.jf -k 4

	This command uses a Jellyfish dictionary containing information about high abundance kmers in the genome of interest to screen probes. (We have provided `sp.jf` as an example for you to test the python script, which should pass all 12 probes into the file `3_probes_18_4.bed` . However, you will need to generate your own Jellyfish dictionary for your desired genome in the real case!) To see additional command line arguments available for this script, you can run the python file with the `-h` argument (i.e. `python kmerFilter.py -h').

5. To convert your probe set to their reverse complements, you can use the `probeRC.py` script:

		python probeRC.py -f 3_probes.bed

	This creates a file, `3_probes_RC.bed` containing the reverse complements of all sequences in the original .bed file. To see additional command line arguments available for this script, you can run the python file with the `-h` argument (i.e. `python probeRC.py -h').

6. [Optiona] You can check for secondary structures of probes by calling NUPACK using the `structureCheck.py` script:

		python structureCheck.py -f 3_probes.bed -t 0.4

	This command should pass 6 of 12 example candidate probes. Additional information can be seen in the produced `3_probes_sC.bed` file. To see additional command line arguments available for this script, you can run the python file with the `-h` argument (i.e. `python probeTm.py -h').

7. [Optional] To generate a list of melting temperatures for a given probe set, you can use the`probeTm.py` script:

		python probeTm.py

	or

		python probeTm.py -f 3.txt

	The first command will allow you to enter a sequence interactively to retrieve its computed melting temperature. The second command takes a two column .txt file with the sequence in column 2 (tab delimited) and outputs a new file (`3_tm.txt`) with a 3rd column of Tms.


That's all! If you made it through these all without any errors thrown about missing dependencies or modules, you are all set to run OligoMiner on your own computer. Happy FISHing!

### Notes on running OligoMiner on new genomes

You'll need to download your genome of interest in FASTA format and prepare index/dictionary files for your NGS aligner and optionally Jellyfish. We recommend using unmasked files for dictionary file construction and repeat-masked files as the input files for `blockParse.py`

## Citation

Please cite according to the enclosed [citation.bib](./citation.bib):

```
@article{Beliveau2018,
        doi = {10.1073/pnas.1714530115},
        url = {https://doi.org/10.1073%2Fpnas.1714530115},
        year = 2018,
        month = {feb},
        publisher = {Proceedings of the National Academy of Sciences},
        volume = {115},
        number = {10},
        pages = {E2183--E2192},
        author = {Brian J. Beliveau and Jocelyn Y. Kishi and Guy Nir and Hiroshi M. Sasaki and Sinem K. Saka and Son C. Nguyen and Chao-ting Wu and Peng Yin},
        title = {{OligoMiner} provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes},
        journal = {Proceedings of the National Academy of Sciences}
}
```

## Questions

Please reach out to [Brian](mailto:beliveau@uw.edu) with any questions about installing and running the scripts, or [open an issue](../../issues/new) on GitHub.

## License

We provide this open source software without any warranty under the [MIT license](https://opensource.org/licenses/MIT).

## Contributing

We welcome commits from researchers who wish to improve our software. Please follow the [git flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model. Make all changes to a topic branch off the branch `dev`. Merge the topic branch into `dev` first (preferably using `--no-ff`) and ensure everything works. Code will _only_ merged into `master` for release builds. Hotfixes should be developed and tested in a separate branch off `master`, and a new release should be generated immediately after the hotfix is merged.
