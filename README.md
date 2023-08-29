# RustMelt5
Rust port of MELTING 5

This code is a direct transfer of the logic in MELTING 5 for calculting the melting temperature between an mRNA / RNA duplex. 

The method implemented is from Turner et al. (2006) Nucleic acids research 34: 3609-3614.

## Building
```cargo build --release``` to build an optimized binary.

## Running
```./rust_melt_5 -s RNA_SEQUENCE -n NA_CONCENTRATION -r RNA_CONCENTRATION```

## Motivation
Some projects currently use `rmelting` which serves as R bindings for the MELTING 5 tool, which contains many methods for calculting oligonucleotide melting temperature. Unfortunately, MELTING 5 is quite slow for various reasons, and the R bindings add even more overhead. Below is an example of a run using *just* the source Java MELTING 5 jar compared to the optimized Rust build.
![test](https://github.com/acrinklaw/RustMelt5/assets/25108033/5f4ca91b-951e-416b-884b-51b8da1f461e)


144X speedup! If we have 3000 oligos that need to be run on a moderately sized gene, rmelting is usually around .3s per oligo, resulting in 15 minutes of real time computation. The Rust version, at a modest run time of .001s, would take 3 seconds to perform equivalent computations.

## Future Plans
Explore adding [PyO3](https://pyo3.rs/v0.19.2/) support to enable us to have Python bindings to the function and avoid having to call an external executable. Overhead of PyO3 would need to be explored.

Implement methods for more hybridization types, beyond just mRNA/RNA

## How It Works
The melting temperature of an mRNA / RNA duplex is fairly straightfoward. There are empirically determined values for the entropy / enthalpy of specific amino acid pairs across the length of the duplex, as well as a factor for the initiation of hybridization. There are a few corrections that are applied specifically for this kind of duplex.

The first correction is a transfer from entropy calculations in 0.1M Na+ to 1M Na+. The second correction is to account for the concentration of Na+ that is actually present in the solution. There is a third correction that is only considered when there is a terminal Adenine or terminal Uracil. (Seems to be thermodynamically unfavorable)

1. Calculate the sum of entropy and enthalpy for amino acid pairs in the sequence. This is using the nearest-neighbor method so we have lookups for pairs such as AA / UU where UU would be on the complementary strand. Since we are assuming perfect matches this makes things easy and we don't truly need to compute the complementary sequence.
2. Apply entropy corrections for 1M Na+ and then true Na+ concentrations in solution.
3. Apply terminal A/U penalty if applicable. 
4. Add initiation values.
5. Compute melting temp: `let tm = enthalpy / (entropy + 1.99 * (nucleotide_concentration/ 4.0).ln()) - 273.15;`
