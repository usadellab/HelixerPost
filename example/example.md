# Example HelixerPost data

## source
The following was created from the first 21384 base pairs
of NC\_023163.2 in the Oryza\_brachyantha genome, version
GCF\_000231095.2\_ObraRS2 from refseq.

The genome was preprocessed into numeric form (genome\_data.h5) with `fasta2h5py`,
and base-wise genic predictions (predictions.h5) were created with `HybridModel.h5`
both from the Helixer project (https://github.com/weberlab-hhu/Helixer).

## usage
From this directory:

```
helixer_post_bin genome_data.h5 predictions.h5 output.gff3
```

The same result is obtained by explicitly passing the bundled config file
(which holds the binary's built-in defaults — see `default_config.yml`):

```
helixer_post_bin --config default_config.yml genome_data.h5 predictions.h5 output.gff3
```

A rating / stats file is also produced. By default it is written alongside the
GFF as `<gff>.rating`; pass `--rating <path>` to override.

CLI flags override values from `--config` which override the built-in defaults.
Run `helixer_post_bin --help` for the full flag list, and
`helixer_post_bin --print-default-config` to dump a starting YAML template.

## expected output

Using the provided data, the above should create the file 'output.gff3',
and the contents there of are expected to be as follows (as of commit b674af9):

```
##gff-version 3.2.1
##species Oryza_brachyantha
# f4f3f7bf629d02431241bfef18cb8ade  /home/ali/.local/share/Helixer/models/land_plant/land_plant_v0.3_a_0100.h5
##sequence-region NC_023163.2:1-21384 1 21384
NC_023163.2:1-21384	Helixer	gene	17026	18898	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001
NC_023163.2:1-21384	Helixer	mRNA	17026	18898	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001
NC_023163.2:1-21384	Helixer	exon	17026	17356	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.exon.1;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
NC_023163.2:1-21384	Helixer	five_prime_UTR	17026	17115	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.five_prime_UTR.1;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
NC_023163.2:1-21384	Helixer	CDS	17116	17356	.	+	0	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.CDS.1;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
NC_023163.2:1-21384	Helixer	exon	17707	18898	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.exon.2;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
NC_023163.2:1-21384	Helixer	CDS	17707	18305	.	+	2	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.CDS.2;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
NC_023163.2:1-21384	Helixer	three_prime_UTR	18306	18898	.	+	.	ID=Oryza_brachyantha_NC_023163.2:1-21384_000001.1.three_prime_UTR.1;Parent=Oryza_brachyantha_NC_023163.2:1-21384_000001.1
```

