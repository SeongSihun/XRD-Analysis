# CTR.py !!

## Usage

```python
from CTR import *
```
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
```
```python
# The following have been defined in CTR.
#    pi = np.pi
#   inf = np.inf
#   vec = lambda *args: np.array(args)
# class : ['Xray', 'Xray2d', 'Atom', 'Molecule', 'SC', 'FCC', 'BCC', 'Perovskite', 'Film', 'Sample']
```
<br>

### Xray Setting
```python
XRD = Xray(wavelength=Xray.CuKa1) # 1.5406 [Å unit]
Q = XRD.Q(1,0,0) # (100) Scan
E = XRD.Energy
TTHETA = XRD.TTHETA
```

### Atoms
```python
# using atomic number Z
Y  = Atom(39)
Al = Atom(13)
W  = Atom(74)
O  = Atom(8)
```

### Molecules
```python
sto = Perovskite(
    abc = 3.905 * vec(1, 1, 1),
    ABO = (Sr, Ti, O)  # Atoms
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
sro_sto = sro/sto(0,0,1) # sro on sto along with c-axis
wo3_yap = wo3/yap(1,0,0) # wo3 on yap along with a-axis
```

### Films
```python
# Bulk
WO3_Bulk = Film(wo3, (inf,1,1))
YAP_Bulk = Film(yap, (inf,1,1))
# Film
WO3 = Film(wo3, (10,1,1))
YAP = Film(yap, (10,1,1))
# Strained Film
WO3_YAP = Film(wo3/yap(1,0,0), (10,1,1))
```

### Samples
```python
# 2 ways
WO3YAP = Sample(WO3_YAP, YAP_Bulk, nref=(1,0,0))
NNONGOSTO = Sample(NNO, NGO, STO_Bulk, nref=(1,0,0))
# New Version
WO3YAP = WO3/YAP(1,0,0)
NNONGOSTO = NNO/NGO/STO_Bulk(1,0,0)
```

<br>

---

### Scan

#### Atom scan (AFF)
```python
plt.plot(TTHETA, np.abs( Y.aff(Q, E)))
plt.plot(TTHETA, np.abs(Al.aff(Q, E)))
plt.plot(TTHETA, np.abs( O.aff(Q, E)))
```

#### Molecule scan (SF)
```python
plt.plot(TTHETA, np.abs(sto.SF(Q, E)))
plt.plot(TTHETA, np.abs(yap.SF(Q, E)))
plt.plot(TTHETA, np.abs(wo3.SF(Q, E)))
```

#### Film scan (Intensity)
```python
# 2 ways
plt.semilogy(TTHETA, WO3.I(Q, E))
plt.semilogy(TTHETA, YAP_Bulk.I(Q, E))
# New Version
plt.semilogy(XRD.TTHETA, XRD.I(WO3))
plt.semilogy(XRD.TTHETA, XRD.I(YAP_Bulk))
```

#### Sample scan (Intensity)
```python
# 2 ways
plt.semilogy(TTHETA, WO3YAP.I(Q, E))
# New Version
plt.semilogy(XRD.TTHETA, XRD.I(WO3YAP))
```
<br>

---

### 2D Scan (In-plane Data)
```python

XRD2D = Xray2d()
XRD2D.Q(1,0,0)  # In-plane(100) == (0kl) scan
X, Y, I = XRD2D.X, XRD2D.Y, XRD2D.I(WO3/YAP_Bulk(1,0,0))
plt.pcolor(X, Y, I, norm=LogNorm())
```