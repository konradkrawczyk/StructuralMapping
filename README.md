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

### A.2 Step 2 out of 2: Perform the structural mapping.

After all the sequences have been numbered, you structurally map the results by typing:

```
python StructuralAlignment.py structurallymap [experiment name]
```

Where [experiment name] is the name you have given to your experiment in the previous step. Following the 'sample' exmple, this would be:

```
python StructuralAlignment.py structurallymap sample
```

This will save the results in data/structuralmap/[experiment name].

### B Results

Results stored under data/structuralmap/[experiment name] are a collection of json-formatted files with structural annotation for each sequences submitted. A given file is a list from sequence identifiers (as you have specified them in input) to results.

For instance a single sequence has the following annotations:

```
"B_d075d1700ee60c141873c2069c4c2fdb": { <<--- ID of the sequence
		"frame": { <<--- Best alignment results for Chothia framework region.
			"best_pdb": "5v7rH", <<--- Best Chothia framework  region template (PDB + chain)
			"best_sid": 87  <<--- Best Chothia framework  sequence identity (for template above)
		},
		"full": { <<--- Best alignment results for entire variable region
			"best_pdb": "5v7rH", <<--- Best full variable region template (PDB + chain)
			"best_sid": 76  <<--- Best Chothia framework  sequence identity (for template above)
		},
		"fread": { <<--- Best results for each specific CDR
			"H2": { <<---CDR identifier
				"qu": "SGDEGY",  <<---query sequence
				"seq": "NGNSGY", <<---best sequence match according to FREAD -- None if none are found.
				"str": "4n9gH", <<---structure the template is coming from.
				"scr": 25 <<---score of this best template -- the higher the better.
			},
			"H3": {
				"qu": "GSNDWYGIDY",
				"seq": "DPLEYYGMDY",
				"str": "1qkzH",
				"scr": 35
			},
			"H1": {
				"qu": "GFTFSNF",
				"seq": "GFTFTNY",
				"str": "1e4xI",
				"scr": 49
			}
		},
```
### C Parallelization of the protocol.

As you might've noticed, running this pipeline on a single thread on millions of sequences might not be optimal.

For this reason for both steps 1 and 2, we have created utilities which create batch files for processing (multiple instances of the commands above which execute on different portions of the datasets).

### C.1 Parallelize Chothia numbering.

To get the parallelization script for step 1, Chothia numbering, the command is:

```
python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs
```

This will print out on the command line, you might want to pipe it into a suitable bash script:

python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs > number_sequences.sh

You need to specify how many sequences there are in the fasta file, and how many parallel jobs you want to run at the same time. For instance, assuming that I have 6,666,666 sequences and want to run 78 parallel jobs on the sample experiment above, I'd type:

```
python DataProcessing.py parallel sample ../data/raw/fasta/sample.fasta  78 6666666 > number_sequences.sh
```

Now simply run the number_sequences.sh script which will kick off the 78 processes:

```
./number_sequences.sh
```

**NB** You might have to make it executable 

```
chmod u+rx number_sequences.sh
```

### C.2 Parallelize Structural mapping .

Logic is similar as with parallel Chothia numbering.

```
python StructuralAlignment.py parallel [experiment name] number-of-cpus
```
As above, you specify the number of CPUs and the name of the experiment. Continuing with our sample experiment this would be (to run on two CPus):

```
python StructuralAlignment.py parallel sample 2 > structurally_map.sh
```

**NB** If number of CPUs is greater than the number of chunks in the data/numbered/[exp name] folder, this will clearly fail. 




