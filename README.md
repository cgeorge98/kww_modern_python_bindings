# kww

This is **NOT** the home page of **libkww**, a small numeric library that provides
the cosine and sine **Fourier transform of Kohlrausch's stretched exponential function**
exp(-t^beta). The original repoistory can be found [here](https://https://jugit.fz-juelich.de/mlz/kww/-/tree/main/).

This Fourier transform is also known as the **Kohlrausch-Williams-Watts function**.
In condensed-matter physics, it is often used to describe frequency-domain
measurements of relaxation processes.

# Math documentation

Joachim Wuttke:\
*Laplace–Fourier Transform of the Stretched Exponential Function:
Analytic Error Bounds, Double Exponential Transform, and Open-Source Implementation "libkww".*\
[Algorithms 5, 604-628 (2012)](http://dx.doi.org/10.3390/a5040604).

# User documentation

**Author:** Joachim Wuttke, Jülich Centre for Neutron Science (2009-12).\
**Citation:** Please cite the Algorithms paper indicated above.\
**Installation:** CMake based. See file INSTALL in the source archive.\
**License:** GNU General Public License (GPL) version 3 or higher.\
**Manual:** [kww(3)](http://apps.jcns.fz-juelich.de/man/kww.html).

## Wrappers for other programming languages

**Python:** See bindings/python/README in the source distribution.
