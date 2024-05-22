# CTR.py !!

## Usage

### Xray Setting
```
XRD = Xray(wavelength=Xray.CuKa1) # 1.5406 [Å unit]
Q = XRD.Q(1,0,0) # (100) Scan
E = XRD.Energy
TTHETA = XRD.TTHETA
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
sto = Perovskite(
    abc = 3.905 * vec(1, 1, 1),
    ABO = (Sr, Ti, O)
)
yap = Perovskite(
    abc = Molecule.pseudocubic(5.18, 5.32, 7.37),
    ABO = (Y, Al, O)
)
wo3 = Molecule(
    abc = vec(7.69, 7.69, 7.76)/2,
    structure = [
            W(0.5, 0.5, 0.5),
            *[O(*rj) for rj in (np.ones([3, 3]) - np.eye(3))/2]
        ]
)
# strained molecules (constant volume)
sro_sto = sro/sto 
nno_sto = nno/sto
```

### Films
```
# Bulk
WO3_Bulk = Film(wo3, (inf,1,1))
YAP_Bulk = Film(yap, (inf,1,1))
# Film
WO3 = Film(wo3, (10,1,1))
YAP = Film(yap, (10,1,1))
# Strained Film
WO3_YAP = Film(wo3/yap, (10,1,1))
```

### Samples
```
WO3YAP = Sample(WO3_YAP, YAP_Bulk)
NNONGOSTO = Sample(NNO,NGO,STO_Bulk)
```

---

### Scan

#### Atom scan (AFF)
```
plt.plot(TTHETA, np.abs(Sr.aff(Q, E)), color='#800020', label='Sr')
plt.plot(TTHETA, np.abs(Ti.aff(Q, E)), color='#008060', label='Ti')
plt.plot(TTHETA, np.abs(O.aff(Q, E)), color='#006080', label='O')
```

#### Molecule scan (SF)
```
plt.plot(TTHETA, np.abs(ngo.SF(Q, E)), color='#800020', label='ngo')
plt.plot(TTHETA, np.abs(nno.SF(Q, E)), color='#008060', label='nno')
plt.plot(TTHETA, np.abs(sto.SF(Q, E)), color='#006080', label='sto')
```

#### Film scan (Intensity)
```
plt.semilogy(TTHETA, WO3.I(Q), color='#006080', label=f"WO3 Bulk", linestyle='dashed')
plt.semilogy(TTHETA, YAP_Bulk.I(Q), color='#008060', label=f"YAP Bulk", linestyle='dashed')
```

#### Sample scan (Intensity)
```
plt.semilogy(TTHETA, GEN_WO3YAP(N).I(Q), color='#800020', label=f"WO3/YAP; N={N}")
```

