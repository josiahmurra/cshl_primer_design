
# OptimusPrymer
![OptimusPrymer](./figures/Optimus_Prime.webp)
## Project for designing primer sequences

This repository contains pythons scripts for identifying PCR primers. Input file should be a plain text file with a single column of NCBI accession number. OptimusPrymer will automatically download genbank files for each accession and align them using MUSCLE5. 

Example usage:
```
python3 OptimusPrymer
```

Contents of plain text file:
```
NM_001008221.1
NM_000699.4
NM_020978.4
NM_001386925.1
NM_001008219.3
```

