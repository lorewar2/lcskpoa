# LCSKPOA

LCSKPOA: Enabling banded semi-global partial order alignments via efficient and accurate backbone generation through extended lcsk++.

## Setting up / Requirements

Rust 1.65.0+ should be installed.
Requires Rust nightly for SIMD

## Usage

Download the repository.

```bash
git clone https://github.com/lorewar2/lcskpoa.git
cd LCSKPOA
```
Activate rust nightly for the package 

```bash
rustup override set nightly
```
Replicate the results in paper:

```bash
cargo run --releae
```

## How to cite

If you are using LCSKPOA in your work, please cite:

[LCSKPOA: Enabling banded semi-global partial order alignments via efficient and accurate backbone generation through extended lcsk++](https://link.springer.com/article/10.1186/s12859-025-06293-z)
