# 안녕하세요 교수님

- 깃허브로 옮겼습니다.


## Usage

### Sample Settings
```
# 원자번호
Sr, Ti, O = (38, 22, 8)
Nd, Ni, _ = (60, 28, 0)
# Sample Objects
STO = Sample(
  structure = Sample.perovskite(A=Sr, B=Ti, O=O),
  abc = np.array([3.905, 3.905, 3.905]) * 1E-10,
  NaNbNc = (10, 10, 20),
  nhkl= np.array([0, 0, 1])
)
NNO = Sample(
  structure = Sample.perovskite(Nd,Ni,O),
  abc = STO.const_volume_abc(Sample.pseudocubic(np.array([5.387, 5.383, 7.610]) * 1E-10)),
  NaNbNc = (10, 10, 1),
  nhkl= np.array([0, 0, 1])
)
```

### Detector Settings
```
# Detector Object
detector = Detector(
  wavelength   = 1.5406e-10,
  samples      = (STO, NNO),
  pixel_length = 0.005,
  size         = (100, 100)
)

```

### Experiment Setting (Condition)
```
(alpha, phi) = 0.1, 0
hkl = np.array([0, 1, 2])
```

### Result
II = detector.scan(hkl, (alpha, phi))
