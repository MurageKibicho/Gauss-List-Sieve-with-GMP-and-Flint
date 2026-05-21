# Lattice Sieve in C
This is a rewrite of [Orry Gooberman's Gaussian Sieve](https://github.com/orrygoob/Gauss-Sieve/blob/main/README.md) for lattice reduction. 

I wrote about the implementation on [LeetArxiv](https://leetarxiv.substack.com/p/gauss-lll-sieve)
Changes:
1. C not C++.
2. Big Integer support using GMP, MPFR and Flint.

lol i hardcoded the matric from [SVP 10 LLL Matrix](https://github.com/orrygoob/Gauss-Sieve/blob/main/challenges/svp10LLL.txt).

Clone the repo and run the code using FLINT 3.4 and the command:

```
clear && gcc GaussSieve.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o
```

The output is a database of short vectors. (It's what I needed).

```
Basis Vectors
-30,-277,-297,269,384,20,409,-328,66,-232,
-86,-313,25,244,-451,-299,289,-31,-163,-533,
496,129,-203,348,-248,-354,-453,-415,-114,93,
-8,-317,-431,-859,-494,-217,-50,231,-105,9,
196,-146,233,604,-57,-53,42,115,839,525,
-185,610,244,-57,-200,329,653,-37,-625,347,
65,183,-253,114,-214,712,-644,490,243,-372,
-1036,-185,-173,-213,144,-543,-143,-52,396,-179,
275,89,447,-666,283,-7,50,-282,632,-559,
348,430,-871,235,-77,-570,453,337,-243,55,

Vectors in Database (Norm, PreImage, Lattice point)

DB(0):0
PreImage: 0,0,0,0,0,0,0,0,0,0,
Lattice point: 0,0,0,0,0,0,0,0,0,0,

DB(1):1831551
PreImage: 0,0,0,0,0,0,0,0,0,1,
Lattice point: 348,430,-871,235,-77,-570,453,337,-243,55,

DB(2):1600978
PreImage: 0,0,0,0,0,0,0,0,1,0,
Lattice point: 275,89,447,-666,283,-7,50,-282,632,-559,

DB(3):1710414
PreImage: 0,0,0,0,0,0,0,1,0,0,
Lattice point: -1036,-185,-173,-213,144,-543,-143,-52,396,-179,

DB(4):1519728
PreImage: 0,0,0,0,0,0,1,0,0,0,
Lattice point: 65,183,-253,114,-214,712,-644,490,243,-372,

DB(5):1556163
PreImage: 0,0,0,0,0,1,0,0,0,0,
Lattice point: -185,610,244,-57,-200,329,653,-37,-625,347,

DB(6):1479430
PreImage: 0,0,0,0,1,0,0,0,0,0,
Lattice point: 196,-146,233,604,-57,-53,42,115,839,525,

DB(7):1382287
PreImage: 0,0,0,1,0,0,0,0,0,0,
Lattice point: -8,-317,-431,-859,-494,-217,-50,231,-105,9,

DB(8):1010869
PreImage: 0,0,1,0,0,0,0,0,0,0,
Lattice point: 496,129,-203,348,-248,-354,-453,-415,-114,93,

DB(9):853468
PreImage: 0,1,0,0,0,0,0,0,0,0,
Lattice point: -86,-313,25,244,-451,-299,289,-31,-163,-533,

DB(10):719100
PreImage: 1,0,0,0,0,0,0,0,0,0,
Lattice point: -30,-277,-297,269,384,20,409,-328,66,-232,

DB(11):1867028
PreImage: -1,0,-1,-1,0,-1,0,0,0,1,
Lattice point: 75,285,-184,534,481,-348,-106,886,535,-162,

DB(12):1786512
PreImage: 1,0,0,1,1,1,0,0,0,0,
Lattice point: -27,-130,-251,-43,-367,79,1054,-19,175,649,

DB(13):1714124
PreImage: 0,0,1,0,0,1,0,1,0,0,
Lattice point: -725,554,-132,78,-304,-568,57,-504,-343,261,

DB(14):2233034
PreImage: 1,0,1,1,1,2,1,1,1,0,
Lattice point: -412,696,-189,-517,-602,216,517,-315,707,-21,
```


All credit goes to [Orry Gooberman's Gauss Sieve](https://github.com/orrygoob/Gauss-Sieve/blob/main/README.md). His code is super articulate.

