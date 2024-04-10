# CTR Analysis (VER 3)
## Available Functions (~ Now)
- Atomic Form Factor
- Structure Factor
<!-- - available on epitaxial film (of which every nhkl have the same orientation as the substrate.)
- available only where nhkl (by pseudo) is parallel to (001), (010), or (100). -->

## Usage

### Xray Setting
```
CuKa1   = Xray(Xray.CuKa1)
# CuKa1 = Xray(wavelength = 1.5406) # Å unit
```

### Atoms
```
# using atomic number Z
Sr = Atom(38)
Ti = Atom(22)
O  = Atom(8)
Ru = Atom(44)
```

### Molecules
```
# Å unit
sto = Molecule(
    structure = Molecule.perovskite(Sr, Ti, O),
    abc = (3.905, 3.905, 3.905)
)
sro = Molecule(
    structure = Molecule.perovskite(Sr, Ru, O),
    abc = Molecule.pseudocubic(5.567, 5.5304, 7.8446)
)
# strained molecules (constant volume)
sro_sto = sro/sto   # = sto.strain(sro)
nno_sto = nno/sto   # = sto.strain(nno)
```

### Films
```
STO = Film(sto, (10, 10, 100))
SRO = Film(sro, (10, 10, 10))
#
SRO_sto = Film(sro/sto, (10, 10, 100))
```

### Samples
```
SRO_STO = Sample(SRO_sto/STO)
NNO_SRO_STO = Sample(NNO/SRO/STO)
```

---

### Detector
```
XRD = Detector(Xray, [Atom], [Molecule], [Film], [Sample])
```

#### AFF scan
```
XRD = Detector(CuKa1, atom = [Sr, Ti, O])
AFF = XRD.AFF(REF=True)
#
plt.plot(XRD.TTHETA, AFF[0], label='Sr')
```

#### SF scan
```
XRD = Detector(CuKa1, molecule=[sro, sro_sto])
XRD.align(nref = [sro(0,0,1), sro_sto(0,0,1)])
SF  = XRD.SF()
#
plt.plot(XRD.TTHETA, np.abs(SF[0]), label='sro', linestyle='dashed')
```
