TABLE2ASN AUTOMATED SEQUENCE SUBMISSION PROGRAM

table2asn is a command-line program that automates the creation of sequence records for submission to GenBank. It uses many of the same functions as Genome Workbench, but is driven entirely by data files, and the records it produces need no additional manual editing before submission.

Entire genomes, consisting of many chromosomes with feature annotation, can be processed using table2asn and properly prepared input files.

In general, table2asn will recognize files with the same basename as the input sequence file, and will output an ASN.1 (Abstract Syntax Notation 1) text file with the same basename and a .sqn suffix. When additional arguments are included additional output files can be generated. For example: validation files with .val suffix and a summary of the .val files output in a file with .stats suffix, a discrepancy report with .dr suffix, a GenBank flatfile with .gbf suffix.

table2asn is supplied for a variety of popular computer platforms, and is available by anonymous ftp at:

  ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn

To get a summary of the command-line arguments that table2asn has at its disposal, run:

  table2asn -help

COMMON SUBMISSION REQUIREMENTS

table2asn expects a file of nucleotide sequence data, in FASTA format, with a .fsa suffix, but it can also extract sequence information from ASN.1 objects. For submitting a single sequence record, the file is selected with the -i argument:

  table2asn -i sde3.fsa

Other types of data files for this record will have the same base name (sde3), but with different suffixes. For example, a feature table file, in a five-column format described later, would have a .tbl suffix. Running table2asn -i sde3.fsa will automatically look for sde3.tbl without needing to explicitly specify it in another argument.

(The .fsa file could also contain a set of concatenated FASTA sequence records, in which case the .tbl file could contain a matching set of concatenated feature tables, as discussed below.)

Every submission also needs a template file containing contact information (for a person to whom questions about the submission can be addressed) and a submission citation (which lists the authors who get scientific credit for the sequencing). This file is in text ASN.1 format, has a .sbt suffix, and, since its name is not likely to correspond to any sequence record, it is selected with the -t argument:

  table2asn -t sangerLab.sbt -i sde3.fsa

In addition, information about the biological source is required. This information can be provided via the definition line of the FASTA file, the command line with the -j argument, or a .src file, as described below.

The output file is in text ASN.1 format, with a .sqn suffix (for Sequin, the predecessor of Genome Workbench) and the same base name, in this case generating sde3.sqn. This output name can be overridden with the -o argument, if necessary:

  table2asn -t sangerLab.sbt -i sde3.fsa -o helicase.sqn

The .sqn file is generally suitable for submission to GenBank. The -binary argument will produce binary ASN.1, if desired.

INPUT FILE TYPES AND SUFFIXES

Details about the various input file types are given in the appropriate sections, but a quick overview is shown here.

Nucleotide sequence data without feature annotation can be obtained from FASTA (.fsa) files.

Biological features for the sequence can be read from Feature Table (.tbl) and GenBank-specific GFF3 (.gff) files.

Sequence and feature data can be read in a single step from ASN.1 (.asn or .sqn) files.

Correction of CDS translation discrepancies (due, for example, to post-transcriptional RNA editing) can be made with sequences in .pep files, finding the appropriate CDS feature by matching the protein_id qualifier.

Source qualifiers for sets of related sequences can be stored in Source Table (.src) files. [Alternatively, source qualifiers can be provided in the definition line of the FASTA file, or in the command line with the -j argument.]

Base calling reliability scores can reside in Quality Score (.qvl) files.

The submission information is usually read from a .sbt template file containing a Submit-block object. It can also be obtained, along with sequence and feature data, from an earlier .sqn file containing a Seq-submit object, in which case a template file indicated by the -t argument would not be needed.

SUBMISSION TEMPLATE WEB SITE

A submission template file can be created and then downloaded by filling in the web form at:

  https://www.ncbi.nlm.nih.gov/WebSub/template.cgi/

Templates can also be exported from the Genome Workbench application.

BULK SUBMISSION VARIANT

For bulk submission of multiple records, all input files are placed in a single directory. The path to this directory is specified in the -indir argument. When -indir is present, each nucleotide file is identified by suffix and processed, one at a time, along with its associated data files.

The nucleotide input file suffix, normally .fsa, can be changed by the -x argument (e.g., to .asn). The output file suffix, normally .sqn, can be changed by the -out-suffix argument.

The -E argument causes recursive exploration of all subdirectories within the input directory.

Output files are placed with the input files, unless overridden with -outdir. In this case, when -M n or -V v or -Z is used, the name of the output directory is the basename of the .stats and .dr files that are generated (those arguments are discussed in more detail below). 

A 5-column feature table or GenBank-specific GFF3 file can be explicitly selected with the -f argument, but not when -indir is used.

Most other command-line customization arguments behave the same whether -i or -indir is used to select sequence records.

ASN.1 FORMAT RECORDS

Sequence files in ASN.1 format may be the product of an earlier submission (using table2asn or Genome Workbench), or can be for previously-published records downloaded from NCBI, either from the Entrez web site:

  https://www.ncbi.nlm.nih.gov/nuccore/U54469

or by using the Entrez Utilities API or its command-line front-end, Entrez Direct:

  efetch -db nuccore -id U54469 -format asn > U54469.asn

ASN.1 files for submission may also be prepared by some commercial sequence analysis software packages.

ANNOTATED FASTA FORMAT

table2asn can read nucleotide sequences in FASTA format. FASTA files start with a definition line and continue with lines of sequence letters. The definition line begins with a right angle bracket (">"), followed immediately by a sequence identifier (SeqID) string, and optionally a space followed by a title. See:

  https://www.ncbi.nlm.nih.gov/genbank/fastaformat/

The FASTA title (definition line) can contain bracketed information about the biological source of the organism. table2asn uses this fielded information when building the submission. A sample definition line is:

  >sde3n [organism=Arabidopsis thaliana] [ecotype=C24] [chromosome=1] RNA helicase SDE3

[Alternatively, the source information in the same bracketed format can be included with a -j argument in the table2asn command line, as seen in some examples at the end of this document. This strategy is useful when all the sequences in the FASTA file(s) have the same source information.]

Other common elements include [topology=circular] and [location=mitochondrion]. RNA viruses would be indicated by [molecule=rna] and [moltype=genomic]. Many other source qualifiers, such as map, strain, clone, isolate, cell-line, and cultivar, can also be used. See:

  https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/

For organisms that are not commonly submitted with table2asn, the genetic code can be supplied with [gcode=x]. For example, [gcode=1] for nuclear genes of many eukaryotes, or [gcode=11] for most prokaryotes. Including the correct genetic code will ensure proper translation of protein coding region features.

The mitochondrial and plastid genetic codes for eukaryotes can differ from their nuclear genetic codes, so those sequences require separate attributes. For example, some green algae use [mgcode=22] and [pgcode=11] for their mitochondrial and plastid genetic codes, respectively.

For the complete set of known genetic codes, see:

  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

Note that the definition line must be a single line, with no return or newline characters. Some word-processors will word-wrap text, either during display or when saving to a file, and care must be taken to avoid unwanted newlines introduced by the editor.

A sample FASTA nucleotide sequence, with structured source fields, is:

  >sde3n [organism=Arabidopsis thaliana] [ecotype=C24] [chromosome=1] RNA helicase SDE3
  TCTTATTTATTTGATTAAATATGCATCTGTCCAATTATAATTATTACGGCTCCAAAGGTCAATTCTTAGA
  ATTACCGTTGTAATTAACAACGATTGGGCTTAACATATACCTCTTAGGCCCATAACATATACAAAACCCA
  AAAGGCCCAAATGTTTGATGCACGTGCTTTACTTGAGCGAATTAATTCCCAATTTCAGCAATTTTAATTT
  TTCTACGACTCTTCAGTCTTCACTCACTCTTTTCATGTTTCTTCTCCTTTGAAGCCTGCCTGCGTTAGTC
  TGGCTTCATTGCTTCTCCATTTCTTGGTGTGATCGAATCAAAGAGTGTAACCCATTTTGCTACTGATTCA
  GTACGTATGATCAATTCTCTCAATTTCAGTTAATCTCATGCTCAATTTCGTTTTCTGTGTTAGGGTTTTG
  GGTTTTTGTTATGCTCTGAGTCTAGTTCACGCTACTCGAATTCCAATACAATCCTCTTAGCGTCATGTTG
  ...

NUCLEOTIDE SEQUENCE SETS

In addition to a single FASTA nucleotide sequence per file, multiple individual sequence records can be packaged in a single file. These submissions may consist of the same sequence region from different strains or organisms, or multiple related sequences, such as the sequences of a genome.

The type of sequence set data is specified by the -a argument, which can take the following values:

  a     Any format, including a single FASTA or ASN.1
  s     Batch set of unrelated sequences

The -a "s" argument is used for a file containing independent FASTA records.

The master genome -M n argument incorporates the -a s functionality.

For .fsa files containing multiple independent FASTA records, the corresponding .tbl or .gff files could contain a set of concatenated feature tables, associated by use of matching sequence identifiers, as discussed below.

FEATURE TABLE FORMAT

table2asn reads features from a five-column, tab-delimited file called a feature table, usually with a .tbl suffix.

The first line of the feature table file should start with a right angle bracket followed immedidately by the the word "Feature"", which is followed by a space, a sequence identifier string (SeqID), and, optionally, a space and a feature table name. The SeqID must match that in the corresponding .fsa file. For example:

  >Feature sde3n

Subsequent lines of the table list the features. Columns are separated by tabs.

The feature table specifies the location and type of biological features. table2asn will process the feature intervals and translate coding region (CDS) features into proteins.

The first and second columns are the start and stop locations of the feature, respectively, the third column is the type of feature (the feature key, e.g., gene, mRNA, CDS), the fourth column is a qualifier name (e.g., "product", and the fifth is a qualifier value (e.g., the name of the protein or gene). A simple example is:

  >Feature sde3n
  240     4084    gene
                          gene       SDE3
  240     1361    mRNA
  1450    1641
  1730    3184
  3275    4084
                          product    RNA helicase SDE3
  579     1361    CDS
  1450    1641
  1730    3184
  3275    3880
                          product    RNA helicase SDE3

Prokaryotic and eukaryotic genome submissions impose additional qualifier requirements for gene, mRNA, and CDS features, which are discussed in detail later.

If a feature contains multiple intervals, each interval is listed on a separate line by its start and stop position. Features that are on the complementary strand are indicated by reversing the interval locations. Locations of partial (incomplete) features are indicated with a '>' or '<' next to the number. For partial features column 1 always has '<' and column 2 always has '>', regardless of the strand on which the feature appears.

A simple example  of a gene that is on the minus strand and is partial at its 3' end is:

  >Feature eno3
  1018    >1      gene
                          gene       ENO3
  1018    >1      CDS
                          product    beta-enolase

Gene features should contain a single interval, and the location should cover the intervals of all other features that are considered to be part of that gene or stages of its expression. Genes can only have multiple intervals when they span the origin of a molecule with circular topology, or in cases of trans-splicing, in which case the trans-splicing exception must also be present.

If the gene feature spans the intervals of the CDS or mRNA features for that gene, there is no need to include gene qualifiers on those features in the table, since they will be picked up by overlap. Use of the overlapping gene can be suppressed by adding a gene qualifier with the value "-". This is important when, for example, a tRNA is encoded within an intron of a housekeeping gene.

Translation exception qualifiers in coding regions are parsed from the same style used in the GenBank flatfile:

                          transl_except       (pos:591..593,aa:Sec)

These are typically due to RNA editing, or the presence of suppressor tRNAs, or to 3'-partial termination codons completed by polyadenylation of the mRNA transcript.

The codon recognized and anticodon position of tRNAs can also be given:

                          codon_recognized    TGG
                          anticodon           (pos:7591..7593,aa:Trp)

In addition to the standard qualifiers seen in GenBank format, several other tokens are used to direct values to specific fields in the ASN.1 data. These include gene_syn or gene_synonym, gene_desc, prot_desc, prot_note, region_name, bond_type, and site_type.

Gene Ontology (GO) terms can be indicated with the following qualifiers:

                          go_component        endoplasmic reticulum|0005783
                          go_process          glycolysis and gluconeogenesis|57|89197757|ACT,TEM
                          go_function         excision repair|93||IPD

The value field is separated by vertical bars '|' into a descriptive string, the GO identifier (leading zeroes are retained), and optionally a PubMed ID (or GO Reference number starting with a leading 0) and one or more GO evidence codes. The evidence code is the fourth token, so include blank fields, as necessary.

The full set of available feature qualifiers is shown at:

  https://www.insdc.org/documents/feature_table.html

BIOLOGICAL EXCEPTIONS

Exceptional biological situations can be annotated in Feature Tables by use of an exception qualifier:

                          exception           ribosomal slippage

The most common legal exception qualifier values are:

  RNA editing
  ribosomal slippage
  trans-splicing

The International Nucleotide Sequence Database collaboration allows "RNA editing" to appear in the /exception qualifier, and displays "ribosomal slippage" and "trans-splicing" in their own qualifiers. Any other valid exception (such as "rearrangement required for product" or "nonconsensus splice site") is shown in a /note qualifier in the flatfile, but must be entered as an exception in the feature table.

GENERAL FEATURE FORMAT

GFF3 files in GenBank-specific format:

  https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/

can be also used to annotate new submissions, in which case the identifier in column 1 must be identical to the SeqID in the FASTA file, just as for .tbl files.

Several arguments are relevant when using GFF3 input files, and are in addition to other arguments for genome submissions and gapped sequence assemblies (discussed later):

Always use -J, which ensures proper processing of protein_id/transcript_id qualifier pairs.

Also use -c w, which runs the WGS cleanup function.

For eukaryotes, add -euk, which will create missing mRNA features.

If locus_tags are not in column 9 of every gene in the GFF3 file, use -locus-tag-prefix followed by the registered prefix.

Be sure to include the product in column 9 of the CDS or exon rows. Otherwise, the protein names will be "hypothetical protein".

SOURCE TABLE FORMAT

For sets of sequences, a source qualifier table can optionally be placed in a tab-delimited file with a .src extension. The first line gives the source qualifier names, separated by tabs. The first column must be the sequence identifier:

  sequence_id    organism                    strain

The remaining lines each give the source qualifiers for one sequence:

  zm57           Zea mays                    A69Y
  sc16           Saccharomyces cerevisiae    S288C
  ...

The complete set of source qualifiers is shown at:

  https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html

PEPTIDE SEQUENCE FORMAT

Peptide sequences are FASTA files with a .pep extension that can substitute for the translated product of a CDS feature. This can compensate for an exceptional situation where CDS translation does not produce the established sequence. The FASTA defline with an angle bracket and sequence identifier is required:

  >sde3p
  MSVSUYKSDDEYSVIADKGEIGFIDYQNDGSSGCYNPFDEGPVVVSVPFPFKKEKPQSVTVGETSFDSFT
  VKNTMDEPVDLWTKIYASNPEDSFTLSILKPPSKDSDLKERQCFYETFTLEDRMLEPGDTLTIWVSCKPK
  ...

and the SeqID must match a CDS feature's protein_id qualifier in the original .tbl file:

                          protein_id          lcl|sde3p

The protein_id needs to explicitly use a 'lcl|' prefix before the SeqID string to indicate a local identifier. A local sequence identifier is assumed when reading FASTA, but a database accession is assumed in the feature table.

An appropriate exception, such as "RNA editing", is expected when a .pep file is used since the protein sequence in the final .sqn file will be different from the conceptual translation of the coding region.

QUALITY SCORES FORMAT

Quality scores can be supplied in .qvl files. These generate Seq-graph data that will be attached to the nucleotide sequence from the .fsa file.

  >chr1
  51 63 70 82 82 82 90 90 90 90 86 86 86 86 90 90 90 90 90 86
  86 86 86 86 86 86 86 90 90 90 90 90 90 86 86 78 78 90 90 86
  ...

These values can be extracted from the output files of the Phrap and Consed programs used to process raw data from automated sequencing machines.

RECORD CLEANUP AND VERIFICATION

The -c argument runs certain cleanup operations, taking any combination of the following letters:

  b    Basic cleanup (default)
  e    Extended cleanup
  f    Fix product names
  s    Add exception to short introns
  w    WGS cleanup (only needed when a .gff file is used)
  d    Correct Collection Dates (assume month first)
  D    Correct Collection Dates (assume day first)
  x    Extends ends of features by one or two nucleotides to abut gaps or sequence ends
  -    Avoid cleanup

The -V argument is used to run the sequence record validator and/or GenBank flatfile generator. It can take any combination of these letters:

  v    Validate record, saving output to files with a .val suffix and a summary of the .val files in a file with .stats suffix
  b    Generate GenBank flatfile format in files with a .gbf suffix.

The -Z argument runs the sequence discrepancy report, which looks for subtle inconsistencies within a set of related records, and outputs a file with the .dr suffix. The -euk argument asserts eukaryotic lineage for the discrepancy report tests.

GENOME SUBMISSIONS

Genome submissions have additional requirements. For records with annotation, a registered locus_tag prefix is needed. This can be obtained by preregistration of BioProject and BioSample identifiers. Otherwise, submitters can register each of these during a genome or SRA submission.

Gene features in a 5-column feature table should have unique locus_tag qualifiers with registered prefix values:

                          locus_tag           HLC_00001

Prokaryotic CDS features should have unique protein_id qualifiers:

                          protein_id          sde3p

Eukaryotic mRNA and CDS feature pairs should have matching unique transcript_id and protein_id qualifiers:

                          protein_id          sde3p
                          transcript_id       sde3m

Note that the locus_tag can be the base of protein_id and transcript_id, if desired. For example:

                          protein_id          HLC_00001.cds
                          transcript_id       HLC_00001.mrna

In contrast, a .gff file does not require protein_id and transcript_id, since they can be created from the gene feature's locus_tag and the Parent information of each feature. If adding them explicitly, use this format for transcript_id attributes:

  transcript_id=gnl|{dbname}|{ID}

and one of these for protein_id attributes:

  protein_id=gnl|{dbname}|{ID}
  protein_id=gnl|{dbname}|{ID}|gb|{accession}

Further details are available in the eukaryotic annotation guidelines.

See the REFERENCES section for URLs of web pages with more details on annotated prokaryotic and eukaryotic genome submissions.

Master genomes are indicated by the -M argument, which can take the following letters:

  n    Normal
  t    TSA (transcriptome shotgun assembly)

The -M n ("Normal") choice is used for genome submissions and combines flags, replacing "-a s -V v -c f", and also invokes FATAL calls when -Z is included for running the discrepancy report. The -J flag should be added separately.

The -j argument will add source qualifiers. These qualifier values override any conflicting values read from a file.

The -f argument will use a particular 5-column feature table or GenBank-specific GFF3 file as annotation input.

The -Y argument will use the text in a file, while the -y argument will use the quoted text in the command line, to create a COMMENT.

GAPPED SEQUENCE ASSEMBLIES

Not all sequencing projects produce complete contiguous data. For sequenced segments on the same chromosome, if their biological order and relative orientation are known, they can be assembled into a gapped sequence, which is implemented in a "delta" sequence. This consist of islands of known sequence separated by gaps of undetermined sequence and estimated or uncertain lengths.

In FASTA format, the gaps are represented by contiguous runs of "N" bases interspersed between the sequenced regions. The interpretation of gaps, in particular the number of consecutive Ns that is recognized as a gap, is controlled by additional arguments:

The -gap-type argument (table2asn defaults to "scaffold" gap-types when this argument is not present) must be one of the following values:

  scaffold
  short-arm
  heterochromatin
  centromere
  telomere
  repeat
  contamination
  contig
  unknown (obsolete)
  fragment
  clone
  other (for future use)

The -gaps-min integer value is the minimum number of consecutive Ns recognized as a gap.

The -gaps-unknown argument supplies the exact number of consecutive Ns recognized as a gap with unknown length.

The -l argument indicates the type of evidence used to assert linkage across assembly gaps. Values can be:

  paired-ends
  align-genus
  align-xgenus
  align-trnscpt
  within-clone
  clone-contig
  map
  strobe
  unspecified
  pcr
  proximity-ligation

The -linkage-evidence-file argument takes a file listing linkage evidence for gaps of different lengths.

EXAMPLES

Single non-genome submission: a particular .fsa file, and only 1 sequence in the .fsa file, and the source information is in the definition line of the .fsa file:

  table2asn -t template.sbt -i x.fsa -V v

Batch non-genome submission: a directory that contains .fsa files, and multiple sequences per file, and the source information is in the definition line of the .fsa files:

  table2asn -t template.sbt -indir path_to_files -a s -V v

Genome submission: a directory that contains multiple .fsa files of a single genome, and one or more sequences per file, and the source information is in the definition line of the .fsa files:

  table2asn -t template.sbt -indir path_to_files -M n -Z

Genome submission for the most common gapped situation (runs of 10 or more Ns represent a gap, and there are no gaps of completely unknown size, and the evidence for linkage across the gaps is "paired-ends"), and the source information is in the definition line of the .fsa files:

  table2asn -t template -indir path_to_files -M n -Z - gaps-min 10 -l paired-ends

GFF3 prokaryotic genome example where the GFF3 file does not include locus_tags, and its basename does not match that of the FASTA file:

  table2asn -M n -J -c w -t template.sbt -gaps-min 10 -l paired-ends -locus-tag-prefix XXXX \
    -j "[organism=Escherichia coli] [strain=abcd] [gcode=11]" -i fasta_file -f gff_file \
    -o output_file.sqn -Z

GFF3 eukaryotic genome example where the GFF3 file does include locus_tags, but its basename does not match that of the FASTA file:

  table2asn -M n -J -c w -euk -t template.sbt -gaps-min 10 -l paired-ends \
    -j "[organism=Loa loa] [isolate=F231]" -i fasta_file -f gff_file \
    -o output_file.sqn -Z

REFERENCES

  FASTA format:

    https://www.ncbi.nlm.nih.gov/genbank/fastaformat/

  Source modifiers:

    https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/

  Genetic codes:

    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 
  Prokaryotic and Eukaryotic Genome Submission Guide:

    https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/

  Prokaryotic Genome Annotation:

    https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/

  Eukaryotic Genome Annotation:

    https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission/

  Using a GenBank-specific GFF file:

    https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/

  Validation errors:

    https://www.ncbi.nlm.nih.gov/genbank/genome_validation/

  Disrepancy report:

    https://www.ncbi.nlm.nih.gov/genbank/asndisc/

  Create a BioProject before submission if a locus_tag prefix is needed:

    https://submit.ncbi.nlm.nih.gov/subs/bioproject/

  Create a BioSample before submission if a locus_tag prefix is needed:

    https://submit.ncbi.nlm.nih.gov/subs/biosample/
    