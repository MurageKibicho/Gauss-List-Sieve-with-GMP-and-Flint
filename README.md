# Lattice Sieve in C
This repo refactors [Orry Gooberman's Gauss Sieve](https://github.com/orrygoob/Gauss-Sieve/blob/main/README.md) Gauss sieve for lattice reduction. 
Changes:
1. Use C not C++.
2. Big Integer support using GMP, MPFR and Flint.

lol i hardcoded the matric from [SVP 10 LLL Matrix](https://github.com/orrygoob/Gauss-Sieve/blob/main/challenges/svp10LLL.txt).

You can run the code using FLINT 3.4 and the command.
All credit goes to [Orry Gooberman's Gauss Sieve](https://github.com/orrygoob/Gauss-Sieve/blob/main/README.md). His code is super articulate.
```
//clear && gcc GaussSieve.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o
```
