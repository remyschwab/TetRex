

## Long, Simple Query

### "Ideal"
Large k, large bin size, long, simple query

- Query: ACGTTTCGATATGCGATAGCTATAGCTAGA
- Query Length: 31
- k = 21
- Theoretical FPR: 1.0E-7
- Bin Size: 43273758033

```shell
## k = 21
./bin/kbioreg benchmark -m na -k 21 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A.T.A.G.C.T.A.T.A.G.C.T.A.G.A."
[0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0]
```

No FPs

--------------------------------

### Kmer Length stressing
- Query: ACGTTTCGATATGCGATAGCTATAGCTAGA
- Query Length: 31
- k = 17
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## k = 17
./bin/kbioreg benchmark -m na -k 17 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A.T.A.G.C.T.A.T.A.G.C.T.A.G.A."
[0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,1]
```

NO False Positives!

---

- Query: ACGTTTCGATATGCGATAGCTATAGCTAGA
- Query Length: 31
- k = 11
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## k = 17
./bin/kbioreg benchmark -m na -k 11 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A.T.A.G.C.T.A.T.A.G.C.T.A.G.A."
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
```

ALL false positives

---

- Query: ACGTTTCGATATGCGATAGCTATAGCTAGA
- Query Length: 31
- k = 13
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## k = 13
./bin/kbioreg benchmark -m na -k 13 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A.T.A.G.C.T.A.T.A.G.C.T.A.G.A."
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0]
```

No False Positives!

---

## Short but simple

- Query: ACGTTTCGATATGCG
- Query Length: 16
- k = 13
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## ACGTTTCGATATGCG
./bin/kbioreg benchmark -m na -k 13 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G."
[0,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0]
```

FP COUNT: 15 | FN COUNT: 0

The probability of a false-positive decreases with the number of kmers. A shorter query has less kmers and results in more FPs

---

- Query: ACGTTTCGATATGCGA
- Query Length: 17
- k = 13
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## ACGTTTCGATATGCGA
./bin/kbioreg benchmark -m na -k 13 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A."
[0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,0,1,0,0,0,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0]
```

FP COUNT: 10 | FN COUNT: 0

---

 Query: ACGTTTCGATATGCGAT
- Query Length: 18
- k = 13
- Theoretical FPR: 1.0E-7
- Bin Size:  43273758033

```shell
## ACGTTTCGATATGCGAT
./bin/kbioreg benchmark -m na -k 13 -s 43273758033 -w 30 -p 20 "AC.G.T.T.T.C.G.A.T.A.T.G.C.G.A.T."
[0,0,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0]
```

FP COUNT: 8 | FN COUNT: 0

---

ACGTTTCGATATGCGATA
FP COUNT: 3 | FN COUNT: 0

---

ACGTTTCGATATGCGATAG
FP COUNT: 1 | FN COUNT: 0
