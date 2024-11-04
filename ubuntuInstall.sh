#!/bin/bash

# Check for Rust
if ! command -v rustc &> /dev/null; then
  echo "Error: Rust is not installed. Please install Rust from https://www.rust-lang.org/tools/install"
  exit 1
fi

# Check for libhdf5-dev
if ! dpkg -l libhdf5-dev &> /dev/null; then
  echo "Installing libhdf5-dev..."
  sudo apt install libhdf5-dev
fi

# Install hdf5-lzf support (Optional)
echo "Do you want to install hdf5-lzf support (y/n)?"
read -r install_lzf

if [[ "$install_lzf" =~ ^([Yy]) ]]; then
  # Download h5py
  wget https://pypi.org/packages/source/h/h5py/h5py-3.2.1.tar.gz

  # Extract and build lzf
  tar -xzvf h5py-3.2.1.tar.gz
  cd h5py-3.2.1/lzf/
  gcc -O2 -fPIC -shared -Ilzf lzf/*.c lzf_filter.c -lhdf5 -L/lib/x86_64-linux-gnu/hdf5/serial -o liblzf_filter.so
  sudo mkdir -p /usr/lib/x86_64-linux-gnu/hdf5/plugins
  sudo cp liblzf_filter.so /usr/lib/x86_64-linux-gnu/hdf5/plugins
  cd ../..
  rm -rf h5py-3.2.1.tar.gz
fi

# Clone HelixerPost repository
git clone https://github.com/TonyBolger/HelixerPost.git

# Build HelixerPost in release mode
cd HelixerPost
cargo build --release

# Move the binary to a directory in your PATH (e.g., /usr/local/bin)
echo "Where would you like to install the binary (e.g., /usr/local/bin)?"
read -r install_dir

if [ ! -d "$install_dir" ]; then
  echo "Creating directory: $install_dir"
  sudo mkdir -p "$install_dir"
fi

sudo mv ./target/release/helixer_post_bin "$install_dir"/helixer_post_bin

echo "HelixerPost has been installed successfully!"
