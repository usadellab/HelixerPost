# HelixerPost
## Dependencies
### Rust
A recent version of Rust (see https://www.rust-lang.org/tools/install)

### hd5 libraries
On fedora & co. you will need `hdf5-devel`, while on ubuntu & co. you will need `libhdf5-dev`.

### hd5 lzf support (Skip if unsure)
If you need the 'lzf' compression support (no longer needed on recent Helixer versions), you will need to download it from h5py (https://pypi.org/project/h5py/) and manually build it as a shared library and install to the hdf5 plugins directory. This will also require a C toolchain - the system provided GCC should be fine. 

`tar -xzvf h5py-3.2.1.tar.gz`

`cd h5py-3.2.1/lzf/`

**Fedora**

`gcc -O2 -fPIC -shared -Ilzf lzf/*.c lzf_filter.c -lhdf5 -o liblzf_filter.so`

`sudo mkdir -p /usr/local/hdf5/lib/plugin`

`sudo cp liblzf_filter.so /usr/local/hdf5/lib/plugin`

**Ubuntu**

`gcc -O2 -fPIC -shared -Ilzf -I/usr/include/hdf5/serial/ lzf/*.c lzf_filter.c -lhdf5 -L/lib/x86_64-linux-gnu/hdf5/serial -o liblzf_filter.so`

`sudo mkdir /usr/lib/x86_64-linux-gnu/hdf5/plugins`

`sudo cp liblzf_filter.so /usr/lib/x86_64-linux-gnu/hdf5/plugins`

## Building HelixerPost

`git clone https://github.com/usadellab/HelixerPost.git`

`cd HelixerPost`

`cargo build --release`

The resulting binary is `./target/release/helixer_post_bin`.

Run `./target/release/helixer_post_bin --help` for the full CLI reference; the
basic invocation is:

```
helixer_post_bin [OPTIONS] <genome.h5> <predictions.h5> <output.gff>
```

In order for Helixer to find this binary, it needs to be on the PATH. The easiest way to achieve this is to copy 
the binary to the bin folder in the virtual environment which you previously created for Helixer 
(e.g. `path_to_Helixer/env/bin` )

## Concept
HelixerPost uses a sliding window assessment to determine regions of the genome which are likely gene containing.
This is then followed by a Hidden Markov Model to convert the base class and coding phase predictions within
that window into one or more gene models, while respecting prior biological knowledge regarding start / stop
codons, RNA splicing etc.  
   
To determine the gene-containing windows, a sliding window of the configured width (e.g. 100bp) are assessed 
for intergenic vs genic (UTR/Coding/Intron) content. The candidate gene containing region starts once the mean 
genic score within the window exceeds the edge threshold, and continues until the mean genic score drops below 
that window. The candidate region is accepted if it also contains at least one window with a genic score above 
the required peak threshold.

## Parameters

```
helixer_post_bin [OPTIONS] <genome.h5> <predictions.h5> <output.gff>
```

Positional arguments:

* `genome.h5` — HDF5 genome used as input to Helixer itself.
* `predictions.h5` — HDF5 base-level predictions emitted by Helixer.
* `output.gff` — destination path for the GFF3 annotation.

Common options (run `--help` for the full list):

* `--rating <path>` — precision/recall stats file. Defaults to `<output.gff>.rating`.
* `-j, --threads <N>` — HMM worker threads. Defaults to `std::thread::available_parallelism()`.
* `-c, --config <path>` — YAML config file. CLI flags below override the file; the
  file overrides built-in defaults. See `example/default_config.yml` for the full schema.
* `--print-default-config` — write the built-in defaults to stdout as YAML and exit.

Pipeline tuning (also configurable via the YAML file):

* `--window-size <bp>` — sliding-window width for genic detection (default `100`).
* `--edge-threshold <f>` — mean genic score required to enter/leave a candidate window (default `0.1`).
* `--peak-threshold <f>` — peak genic score required to accept a candidate window (default `0.8`).
* `--min-coding-length <bp>` — drop genes with total CDS shorter than this (default `60`).
* `--hmm-prob-floor <f>` — minimum probability used before negative-log conversion (default `1e-9`).
* `--hmm-phase-retain <f>` — weight given to the predicted phase signal (default `0.20`).
* `--hmm-start-weight <f>` / `--hmm-stop-weight <f>` — multiplicative weights on the start / stop codon base-match penalties (default `1000`).
* `--hmm-donor-weight <f>` / `--hmm-acceptor-weight <f>` — multiplicative weights on splice donor / acceptor base-match penalties (default `1.0`).

The HMM splice-flag booleans and per-donor fixed penalties are config-file-only
(they appear in `example/default_config.yml` and `--print-default-config`).

## Example Usage

Using the defaults shown above:

```
./target/release/helixer_post_bin example/genome_data.h5 example/predictions.h5 example/output.gff
```

Or, equivalently, with all tuning values loaded from a YAML file:

```
./target/release/helixer_post_bin --config example/default_config.yml \
    example/genome_data.h5 example/predictions.h5 example/output.gff
```

See `example/example.md` for the expected output.


