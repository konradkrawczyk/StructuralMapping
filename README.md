# Structural mapping of Next Generation Sequencing Antibody (Ig-seq) repertoires.

This code is a protocol/pipeline for performing large-scale structural annotations of antibody sequencing repertoires. In other words, giving you majority of info that structural modeling would do, without performing structural modeling per se.

## Requirements.

This code is designed to run ona Unix system (Linux or at the very least Os X). You **might** be able to run it on Windows, but I never tested it in such setup as I could not see a reasonable application for it. 

This protocol assumes that the input is a FASTA-formatted file of amino-acid sequences. If you do not know how to translate your input into amino-acid FASTA, please have a look at IGBlast or IMGT V-Quest.

One thing you need to install yourself is Anarci (free software): http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php

Otherwise it is supposed to work out of the box as described below. If it doesn't drop an email to konrad.krawczyk@stats.ox.ac.uk

## Install

Simply git clone the repo:

```
git clone https://github.com/panzerschnitzel/StructuralMapping.git
```

Navigate to the folder that you just cloned. Throughout I am assuming that you are running from the 'code' directory in this repository.

### A Structurally Mapping a single sequence.

Navigate to the 'code' directory and from here, type

```
python StructuralAlignment.py process_single EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD
```

Results are displayed on output.

### B Results


For instance a single sequence has the following annotations:

```
"B_d075d1700ee60c141873c2069c4c2fdb": { <<--- ID of the sequence
		"frame": { <<--- Best Chothia alignment results for Chothia framework region.
			"best_pdb": "5v7rH", <<--- Best Chothia framework  region template (PDB + chain)
			"best_sid": 87  <<--- Best Chothia framework  sequence identity (for template above)
		},
		"full": { <<--- Best Chothia alignment results for entire variable region
			"best_pdb": "5v7rH", <<--- Best full variable region template (PDB + chain)
			"best_sid": 76  <<--- Best Chothia framework  sequence identity (for template above)
		},
		"cdr": { <<--- Best Chothia alignment for all CDRs taken together.
			"best_pdb": "4nikB", <<--- Template which best matches the three CDRs together.
			"best_sid": 56 <<--- Sequence identity of the entire sequence region.
		}
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

We have a sample set located in data/raw/fasta/sample.fasta which we will use as an example for this protocol.  

To get the parallelization script the command is:

```
python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs
```

This will print out on the command line, you might want to pipe it into a suitable bash script:

python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs > structurally_map.sh

You need to specify how many sequences there are in the fasta file, and how many parallel jobs you want to run at the same time. For instance, assuming that I have 6,666,666 sequences and want to run 78 parallel jobs on the sample experiment above, I'd type:

```
python DataProcessing.py parallel sample ../data/raw/fasta/sample.fasta  78 6666666 > structurally_map.sh
```

You might want to look at the contents of structurally_map.sh if you want to figure out how to run the protocol only on chunks of your data -- it is essentially a collection of python scripts that runs on subsets of lines from the fasta file you specify.

Now simply run the structurally_map.sh script which will kick off the 78 processes:

```
./structuarlly_map.sh
```

**NB** You might have to make it executable 

```
chmod u+rx structurally_map.sh
```

Results stored under data/structuralmap/[experiment name] and are a collection of json-formatted files with structural annotation for each sequences submitted (as described in B). A given file is a list from sequence identifiers (as you have specified them in input) to results.


```
As above, you specify the number of CPUs and the name of the experiment. Continuing with our sample experiment this would be (to run on two CPus):

```
python StructuralAlignment.py parallel sample 2 > structurally_map.sh
```

**NB** If number of CPUs is greater than the number of chunks in the data/numbered/[exp name] folder, this will clearly fail. 




