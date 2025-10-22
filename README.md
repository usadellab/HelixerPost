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

The resulting binary is './targets/release/helixer_post_bin'

Running `./target/release/helixer_post_bin` should show the command line parameters, giving:

`HelixerPost <genome.h5> <predictions.h5> <window_size> <edge_thresh> <peak_thresh> <min_coding_length> <output.gff>`

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

`HelixerPost <genome.h5> <predictions.h5> <window_size> <edge_thresh> <peak_thresh> <min_coding_length> <output.gff>`

genome.h5: The path of the HDF5 formatted genome which was used as input to Helixer itself 

predictions.h5: The path of the HDF5 formatted output from Helixer, containing the base-level predictions

window_size: This determines the number of bases averaged during the sliding window approach (e.g. 100bp)

edge_thresh: This threshold specifies the genic score which defines the start / end boundaries of each 
candidate region (e.g. 0.1)

peak_thresh: This threshold specifies the minimum peak genic score required to accept the candidate region 
(e.g. 0.8)

min_coding_length: The output of the HMM is filtered to remove genes with a total coding length shorter than 
this value (e.g. 60)

output.gff: The path for the output GFF file 

## Example Usage

Using 100bp sliding window for genic detection, 0.1 edge threshold, 0.8 peak threshold (before HMM)

60bp minimum CDS length per gene (after HMM)

From the HelixerPost dir:

`./target/release/helixer_post_bin genome_data.h5 predictions.h5 100 0.1 0.8 60 output.gff`


