# Mod-Minimizer-Minimap2

This is a fork of [minimap2](https://github.com/lh3/minimap2) that replaces the minimizer algorithm with the mod-minimizer scheme, based on the paper:

> **The mod-minimizer: a simple and efficient sampling algorithm for long k-mers**  
> Ragnar Groot Koerkamp, Giulio Ermanno Pibiri  
> bioRxiv 2024.05.25.595898; doi: [10.1101/2024.05.25.595898](https://doi.org/10.1101/2024.05.25.595898)

## Changes

This fork modifies the minimap2 code to implement the mod-minimizer algorithm for finding (w,k)-minimizers on DNA sequences. The mod-minimizer scheme provides a simple and efficient sampling algorithm for long k-mers, improving performance in certain applications.

**Note:** This implementation reduces the minimizer density. If you wish to achieve approximately the same minimizer density as Minimap2's default settings, it is recommended to use the flag `-w 8`.

## Mod-Minimizer Algorithm Overview

The mod-minimizer algorithm finds (w,k)-minimizers on a DNA sequence using the following procedure:

**Notation:**

- **tmer**: A newly constructed t-mer by removing the first base and appending the new base.
- **W**: Set of t-mers in the current window.
- **tmer_i**: The i-th t-mer of the window.
- **kmer_i**: The i-th k-mer of the window.
- **h(a)**: Hash of sequence **a**.
- **rc(a)**: Reverse complement of sequence **a**.
- **pos_W(a)**: Position of **a** within window **W** (0-indexed).
- **M**: List of minimizers.

**Procedure:**

Using a sliding window, construct the entering t-mer:

1. **Update Window**: Remove the oldest t-mer from **W**:  
   `W = W \ {W_0}`

2. **Compute t-mer Info**:  
   `info = min( h(tmer), h(rc(tmer)) )`

3. **Add t-mer to Window**:  
   `W = W ∪ {info}`

4. **Find Minimal t-mer**:  
   `min = min(W)`

5. **Compute Position**:  
   `p = pos_W(min) mod w`

6. **Select k-mer**:  
   `kmer = min( h(kmer_p), h(rc(kmer_p)) )`

7. **Update Minimizers**:  
   `M = M ∪ {kmer}`

For more details on the algorithm and its implementation, please refer to the [original paper](https://doi.org/10.1101/2024.05.25.595898).

## Additional Information

For all other details, usage instructions, and documentation, please refer to the original [minimap2 repository](https://github.com/lh3/minimap2).

