# AGouTI - Annotation of Genomic and Transcriptomic Intervals

## Introduction
High-throughput sequencing techniques have become very popular in molecular biology research. In many cases, obtained results are described by positions corresponding to the transcript, gene, or genome. Usually, to infer the biological function of such regions, it is necessary to annotate these regions with overlapping known genomic features, such as genes, transcripts, exons, UTRs, CDSs, etc. AGouTI is a tool designed to annotate any genomic or transcriptomic coordinates using known genome annotation data in GTF or GFF files. 

#### Main features
1. AGouTI works with coordinates describing positions within the genome and within the transcripts,
2. Ability to assign intragenic regions from provided GTF/GFF annotation (UTRs, CDS, etc.) or de novo (5’ part, middle, 3’ part, whole), 
3. Annotation of intervals in standard BED or custom column-based text files (TSV, CSV, etc.) in any non-standard format,  
4. Flexible handling of multiple annotations for a single region, 
5. Flexible selection of GTF/GFF attributes to include in the annotation. 

<hr>
<br>
<br>
<br>

## Installation

Python >= 3.7 is required to run this software.

You can install AGouTI using `pip` as follows (recommended)

`pip install AGouTI`

or by 

`python setup.py install`

<br>

#### Having troubles with an older Python version?

You can easily manage python versions using `conda` and the concept of virtual environments.

`Anaconda` can be downloaded from https://www.anaconda.com/distribution/#download-section (Python 3.x version)

To create a virtual environment with a specified Python version, you can type

`
conda create --name your_virtualenv_name python=3.8
`

Afterward, you need to activate your virtual environment.

`
conda activate your_virtualenv_name
`

Now you can install AGouTI using `pip` or `conda` and perform your analysis.

After your job is finished, you can leave your virtual environment.

`
conda deactivate
`

As an alternative to conda you can use Python's `venv`

<hr>
<br>
<br>
<br>

## Run AGouTI

You can now access <b>AGouTI</b> as follows:

`
agouti --help
`
<br/>

Running AGouTI is a two-step process. First, you need to create a dedicated database based on your annotation file (GTF/GFF3) using the `create_db` module and then annotate your intervals with the features of interest using `annotate`.

`
agouti create_db --help
`

<br/>

`
agouti annotate --help
`

<br>
<br>

### Step 1. Create the database

We have decided to rely on the SQLite database to efficiently store and access annotation data. Thus, the first step of our annotation pipeline is to create such a database using information included in GTF or GFF3 files. All feature names and attributes from the GTF/GFF3 file are automatically converted to lowercase to uniform the feature selection during the annotation step. By default, the initial database is created in RAM. Then it is inspected to provide the user with a list of features and attributes available in the GTF/GFF file and visualize the hierarchy of those features in a graph-based manner. The inspection process is efficient while the initial database is stored in memory (by default), but for low-memory machines, it can also be stored on the hard drive (an SSD drive would be recommended for speed). Memory consumption depends on the number of features and attributes that must be stored. For example, estimated memory use for Gencode annotation files of the human genome is up to 5GB. Finally, the database is written on a hard drive and can be used for annotation. All those steps are done automatically using create_db mode. Example invocations: 

`
agouti create_db -f GTF -a gtf_of_your_choice.gtf -d database_name
`

or 

`
agouti create_db -f GFF3 -a gff3_of_your_choice.gff3 -d database_name
`

##### Required options

<b>-a</b>, <b>--annotation</b> : Input file containing gene annotations in GTF or GFF3 format.

<b>-f</b>, <b>--format</b> : Input file format (GTF or GFF3), can be gzip-compressed.

<b>-d</b>, <b>--db</b> : Name for the output database.

##### Additional options 

<b>-l</b>, <b>--low-ram</b> : Creates the database as a sqlite3 file directly on your disk. By default, the initial database is created in RAM to quickly inspect contents and relations between features and afterward saved on a local drive. Using this option may significantly slow down database creation. Use only when your RAM size is limited in comparison to the expected size of your database.

<b>-i</b>, <b>--infer_genes</b> : Infer genes. Use only with GTF files that do not contain separate lines describing genes. This step might be very time-consuming.

<b>-j</b>, <b>--infer_transcripts</b> :  Infer transcripts. Use only with GTF files that do not contain separate lines describing transcripts. This step might be very time-consuming. 

##### Output
The output files include: 

1. <b>database_name</b> - the SQLite database file
2. <b>database_name.relations</b> -  text file storing relations between feature types. This file is required for annotation with AGouTI. Therefore it must be stored in the same directory as the database file.
3. <b>database_name.attributes_and_features.pickle</b> - python dictionary stored as a pickle file. This file is required for annotation with AGouTI. Therefore it must be stored in the same directory as the database file. 
4. <b>database_name.database.structure.txt</b> - additional text file listing all the features and attributes present in the GTF/GFF3 file and showing relations between them in a tree-like structure


The output consists of a tree-based structure representing the hierarchy of the features in the GTF or GFF3 file and a list of available attributes for each feature type. That information is by default displayed on stdout. 
Having insights into the database contents, users can choose only a subset of features and attributes to annotate the dataset of interest. 

<br>
<br>

### Step 2. Annotate your file
After creating the database, you can run annotation with AGouTI using the `agouti` annotate command. For example, AGouTI can annotate intervals stored in standard BED or any column-based text file containing information about genomic or transcriptomic coordinates (see the `--custom` option). 

The transcriptomic mode enables the annotation of intervals encoded as positions within the transcripts instead of chromosomes. However, it is instrumental in annotating results from transcript-focused analyses (e.g., RBPs binding sites, identification of RNA structural motifs, etc.), it is essential to use the same source and version of annotation files as during the generation of the intervals submitted for annotation. The transcript IDs often change between the annotation releases as the transcript layout is improved (e.g., ENST00000613119.1 -> ENST00000613119.2 -> ENST00000613119.3, etc.). Since the transcript IDs are part of the coordinate system in transcriptomic mode, any difference will result in an annotation error. 

Basic command:

`agouti annotate -i input.bed -d database_name`

##### Required Options

<b>-i</b>, <b>--input</b> : Input file in BED or another column-based format (see --custom). 

<b>-d</b>, <b>--database</b> : Database created by agouti create_db.

##### Additional options

<b>-m</b>, <b>--custom</b> : Specify that the input text file is in custom format, besides BED. It should contain columns with information about feature id (id), chromosome (chr), start (s), and end (e) coordinates. Users can optionally specify a column with strand information (strand); otherwise, AGouTI will set it to '.'. Format should be specified as column indexes (starting from 1), in the following order: "id,chr,s,e,strand" or "id,chr,s,e", e.g. --custom 1,2,4,5,6. The field separator used in your file can be provided using the --separator option.

<b>-p</b>, <b>--separator</b> : Field separator for the --custom option. Default is tabulator.

<b>-b</b>, <b>--first_base_num</b>' : The first base number in the input file (BED/CUSTOM). Either 0 (0-based coordinates) or 1 (1-based coordinates). Default is 0 (standard for genomic BED files).

<b>-n</b>, <b>--header_lines</b> : The number of header lines. 0 by default. If a single header line is present, set this parameter to 1, etc. 

<b>-t</b>, <b>--transcriptomic</b> : Transcriptomic annotation mode. In this mode, transcript IDs from the GTF/GFF3 are expected to be placed in the first column of provided BED file instead of chromosome names. Coordinates in this mode are assumed to reflect positions within the transcript. Optional.

<b>-f</b>, <b>--select_features</b> : Comma-separated list of feature names to be reported, e.g., "mRNA,CDS". Refer to [db_name].database.structure.txt file written during the database creation for a list of valid features. By default, all features are reported. 

<b>-a</b>, <b>--select_attributes</b> : Comma-separated list of attribute names to be reported, e.g., "ID,description". Refer to [db_name].database.structure.txt file written during the database creation for a list of valid attributes. By default, all attributes are reported.

<b>-c</b>, <b>--combine</b> : List of specific feature-attribute combinations to be reported. Desired combination should be specified in the format: 'feature1-attribute1:attribute2,feature2-attribute1', e.g., "mRNA:description,biotype"  for each mRNA, will provide annotation of mRNA description and mRNA biotype.

<b>-s</b>, <b>--strand_specific</b> : Strand-specific search. 

<b>-w</b>, <b>--completly_within</b> : The feature must lie entirely within the GTF/GFF3 feature to be reported. 

<b>-l</b>, <b>--level</b> : Group results on a specific level (e.g., 1 for gene level, 2 for mRNA, tRNA). Must be one of [1,2]. An annotation may be done on gene or transcript levels (level 1 or 2) so that each output line will correspond to the gene or transcript. Please note that --level 1 cannot be combined with --transcriptomic. Default is 2.

<b>-r</b>, <b>--annotate_feature_region</b> : Report region within the GTF or GFF3 feature, which overlaps with an entry from the input file. Designed to work with `--transcriptomic`. Possible values: 

1. `5 prime` -  when the annotated feature starts within the first quarter of gene or transcript and ends in the first half
2. `middle` -  starts and ends within the second and third quarter
3. `3 prime` - starts within the third quarter and ends in the last one
4. `whole` - starts within the first quarter and ends within the last one. The length of the annotated feature does <b>not</b> exceed 90% of transcript or gene length.
5. `full` - starts within the first quarter and ends within the last one. The length of the annotated feature does exceed 90% of transcript or gene length.

<b>--statistics</b>: Calculate additional statistics. Those will be displayed at the end of the software\'s output (starting with #).

<b>--stats_only</b>: Display statistics only.

##### Output
Output is by default displayed on stdout in the form of a self-explanatory `.tsv` table.
<hr>
<br>
<br>
<br>

## Test case 1

All files used in this case are available at https://github.com/zywicki-lab/agouti (`agouti/agouti_pkg/sample_data.tar.gz`).

### Annotate human single nucleotide polymorphisms (SNPs) stored in the BED file (`common_snp.bed`). 

We've downloaded BED file from the UCSC Table Browser (<i>https://genome.ucsc.edu/cgi-bin/hgTables</i>) using the following filters:

<b>Clade:</b> Mammal
<b>Genome:</b> Human
<b>Assembly:</b> Dec. 2013 (GRCh38/hg38)
<b>Group:</b> Variation
<b>Track:</b> Common SNPs(151)
<b>Table:</b> snp151Common
<b>Position:</b> chrX
<b>Output format:</b> BED - browser extensible data

and saved only the first 1000 SNPs in `common_snp.bed`.

Furthermore, we've downloaded gene annotations (only the X chromosome, genome version <i>GRCh38.p13</i>) in the GFF3 file format from the Ensembl database (<i>https://www.ensembl.org</i>) and subtracted to contain only the first 18178 lines.

To run AGouTI, the User needs to create the database based on the contents of the GFF3 file. It can be done using the command

`agouti create_db -f GFF3 -a Homo_sapiens.GRCh38.105.chromosome.X.gff3.gz -d Homo_sapiens.GRCh38.105.chromosome.X.db`

After the job completes, the User can explore the structure and contents of the database and GFF3 file by examining the `Homo_sapiens.GRCh38.105.chromosome.X.db.database.structure.txt`.

Let’s say that we are not interested in pseudogenes, unconfirmed_transcripts, and non-coding RNAs. To annotate SNPs using our database and calculate additional statistics (`--statistics`), we can type

`agouti annotate -d Homo_sapiens.GRCh38.105.chromosome.X.db -i common_snp.bed -f gene,lnc_rna,exon,mrna,five_prime_utr,three_prime_utr,cds --statistics > annotated_snp.tsv`

You can examine the results stored in the `annotated_snp.tsv`

Please note! Some SNPs lie in the intergenic regions, and the closest gene upstream/downstream is marked as `None`, because we operate near the chromosome border, and no gene is annotated upstream or downstream from such SNPs.

To discard SNPs located in the intergenic regions, we can use `grep`, `awk`, or similar tools, e.g.

`grep -v intergenic annotated_snp.tsv > annotated_snp_intergenic_discarded.tsv`

You can also make annotations on gene instead of transcript level (option `-l`)

`agouti annotate -d Homo_sapiens.GRCh38.105.chromosome.X.db -i common_snp.bed -f gene,lnc_rna,exon,mrna,five_prime_utr,three_prime_utr,cds -l 1 --statistics | grep -v intergenic > annotated_snp_l1.tsv`

You can explore the differences by examining the `annotated_snp_l1.tsv` file.

<hr>
<br>
<br>

## Test case 2

All files used in this case are available at https://github.com/zywicki-lab/agouti (`agouti/agouti_pkg/sample_data.tar.gz`).

### Annotate sample results obtained with the missRNA software (`missRNA.tsv`). Records represent coordinates of the small RNA molecules excised from the longer RNA molecules.

We've downloaded human gene annotations (only the 1st chromosome, genome version <i>GRCh38.p13</i>) in the GFF3 file format from the Ensembl database (<i>https://www.ensembl.org</i>) and subtracted to contain only the first 68571 lines.

To create the databasae, run:

`agouti create_db -f GFF3 -a Homo_sapiens.GRCh38.105.chromosome.1.gff3.gz -d Homo_sapiens.GRCh38.105.chr1.db`

To annotate the file in format other than standard BED using transcriptomic mode, use:

`agouti annotate -i missRNA.tsv -d Homo_sapiens.GRCh38.105.chr1.db -t -r -m 2,1,4,5,6 > missRNA-results.tsv`

The self-explanatory output file is stored as `missRNA-results.tsv` in the same directory as sample datasets.

<hr>
<br>
<br>
<br>


### Contribute

If you notice any errors and mistakes or would like to suggest some new features, please use Github's issue tracking system to report them. You are also welcome to send a pull request with your corrections and suggestions.

### License

This project is licensed under the GNU General Public License v3.0 license terms.

Anytree (Copyright (c) 2016 c0fec0de) and gffutils (Copyright (c) 2013 Ryan Dale) packages are distributed with this software to ensure full compatibility. Documentation, authors, license and additional information can be found at https://anytree.readthedocs.io/en/latest/ and https://github.com/daler/gffutils.
