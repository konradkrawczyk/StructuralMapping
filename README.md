# Structural mapping of Next Generation Sequencing Antibody (Ig-seq) repertoires.

This code is a protocol/pipeline for performing large-scale structural annotations of antibody sequencing repertoires.


## Requirements.

This code is designed to run ona Unix system (Linux or at the very least Os X). You **might** be able to run it on Windows, but I never tested it in such setup as I could not see a reasonable application for it.

The code comes with pre-compiled versions of hmmer so that ANARCI can run.

Otherwise it is supposed to work out of the box as described below. If it doesn't drop an email to konrad.krawczyk@stats.ox.ac.uk

## Install

Simply clone the repo:


```
git clone https://github.com/panzerschnitzel/StructuralMapping.git
```

Throughout I am assuming that you are running

##A. Basic Protocol

We will walk through an

###A.1 Step 1 out of 2: Perform Chothia numbering of the dataset.

Navigate to the 'code' directory and from here, numbering the dataset can be achieved by:

```
python DataProcessing.py number_dataset [experiment_name] [fasta location]
```

We have a sample set located in data/raw/fasta/sample.fasta which we will use as an example for this protocol. This protocol assumes that the input is a FASTA-formatted file of amino-acid sequences. If you do not know how to translate your input into amino-acid FASTA, please have a look at IGBlast or IMGT V-Quest. To number this dataset, type:

```
python DataProcessing.py number_dataset sample ../data/raw/fasta/sample.fasta 
```

Numbered results will be stored under data/numbered/[experiment name], so in this case data/numbered/sample. The results are python-serialized files which you can read by using pickle.

**NB**: if you are re-running numbering do clear the data/numbered/[experiment name] directory beforehand, otherwise you will duplicate data.

###A.2 Step 2 out of 2: Perform the structural mapping.

After all the sequences have been numbered, you structurally map the results by typing:

```
python StructuralAlignment.py structurallymap [experiment name]
```

Where [experiment name] is the name you have given to your experiment in the previous step. Following the 'sample' exmple, this would be:

```
python StructuralAlignment.py structurallymap sample
```

This will save the results in data/structuralmap/[experiment name].

###A.3 Results

###B Parallelization of the protocol.

As you might've noticed, running this pipeline on a single thread on millions of sequences might not be optimal.

For this reason for both steps 1 and 2, we have created utilities which create batch files for processing (multiple instances of the commands above which execute on different portions of the datasets).

###B.1 Parallelize Chothia numbering.

To get the parallelization script for step 1, Chothia numbering, the command is:

```
python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs
```

This will print out on the command line, you might want to pipe it into a suitable bash script:

python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs > number_sequences.sh

You need to specify how many sequences there are in the fasta file, and how many parallel jobs you want to run at the same time. For instance, assuming that I have 6,666,666 sequences and want to run 78 parallel jobs on the sample experiment above, I'd type:

```
python DataProcessing.py parallel sample ../data/raw/fasta/sample.fasta  78 6666666
```





