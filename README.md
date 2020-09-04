# ppr

"Parallel PageRank-Nibble" implementation in C / OpenMP.

## Parallelization

Currently, only the outer loop (over seed nodes) is parallelized, but the two inner loops can be parallelized as well.

## Usage

```
mkdir bin
make
./bin/ppr
```

## Changing datasets

`pack-mm.py` is a Python script that converts a Matrix Market graph to binary, which is read by `ppr`.

## References

[Parallel Local Graph Clustering](https://arxiv.org/pdf/1604.07515.pdf) / Shun et al
