# free-sdm

This project is part of the FRB-CESAB working group FREE and contains data and R scripts to model > 12,000 birds (LifeBirds) and mammals (IUCN) species
at a global scale (50km x 50km resolution) under five GCMs coupled with the RCP8.5 under two times slices.


## Usage and workflow

Clone the repository and edit `coreProgram.R` (essentially `path`, `user` and `n_cores` variables). And then, on a shell:

```bash
Rscript ~/free-sdm/coreProgram.R
```

## Authors

- Nicolas Loiseau
- Maya Guégen
- Wilfried Thuiller
- Nicolas Casajus
