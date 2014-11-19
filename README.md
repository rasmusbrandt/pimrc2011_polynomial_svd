pimrc2011_polynomial_svd
====

**pimrc2011_polynomial_svd** is the simulation environment for
> [R. Brandt][rabr5411] and [M. Bengtsson][matben], "[Wideband MIMO Channel Diagonalization in the Time Domain][pimrc2011_diva]", in _Proc. IEEE Int. Symp. Indoor, Mobile Radio Commun. (PIMRC'11)_, 2011, pp. 1958-1962.

Note that this code was written a couple of years ago, for now defunct versions of Matlab. This code comes with ABSOLUTELY NO WARRANTY.

## Abstract
Methods for spatially diagonalizing wideband multiple-input multiple-output
channels using linear finite impulse response (FIR) filters are investigated.
The PSVD approach by applying the PQRD-BC algorithm for approximate singular
value decomposition (SVD) of polynomial matrices is compared to the approach of
performing a set of conventional SVDs in the Discrete Fourier Transform (DFT)
domain, in terms of complexity and approximation error. Reduced order filters,
based on the DFT-SVDs, are then obtained by optimizing the phases of the
filters. Applying the phase optimized filters as linear filters then forms a
benchmark on the accuracy attainable for any PSVD factorization, for the given
filter length.

Simulations show that the DFT-SVD method has significantly lower complexity than
the PSVD by PQRD-BC, but results in higher order filters. On the other hand,
the PSVD by PQRD-BC yields filters which are close to being perfectly unitary
for all frequencies. To achieve good performance, the reduced order filters are
around one order of magnitude longer than the channel impulse response length.
Therefore there is no gain in performing time domain diagonalization using a
polynomial SVD, compared to using a multicarrier solution.

## Simulation code features

* [Matlab][matlab] implementation of PSVD and PQRD algorithms as proposed in
  [[2]][psvd] and as used in [[1]][pimrc2011_diva] and [[3]][thesis_diva].
* `pmatlib` with `PMatrix` Matlab class, for easy manipulation of polynomial
  matrices in Matlab.
* Matlab implementation of the phase optimization technique proposed in
  [[4]][palomar].

## Running the simulations

1. Ensure that you are able to build MEX files in Matlab. See guide
[here](http://se.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html).
2. For improved performance, build the MEX files:
  * `cd pmatlib` and then `convmat_build`
     The MEX implementation of `convmat` is a direct convolution of the two
     polynomial matrices. If the MEX implementation is not available, for
     example due to not being compiled, a frequency domain Matlab implementation
     is used for the convolution.
  * There is also a MEX file implementing the subproblems of [[4]][palomar]. This
    MEX file however segfaults on the author's new computer, at the time of this
    writing (Nov. 2014). If you are feeling lucky, you may try using this MEX
    file by `cd palomar` and then `palomar_am_loop_build`. If this crashes your
    Matlab, just remove the `palomar_am_loop` MEX file.
3. `matlabpool X` where `X` is the number of cores available
4. Run any `_run` scripts to run the simulation, and the corresponding 
   `_plot` to plot the results.

## pmatlib
As part of this package is the `pmatlib` library, which is a [Matlab][matlab]
library for handling polynomial matrices. By using the provided Matlab class
`PMatrix`, it becomes easy to apply polynomial matrix operations. The `pmatlib`
also contains implementations of the polynomial matrix algorithms in
[[1]][pimrc2011_diva], [[2]][psvd], and [[3]][thesis_diva].

## License and referencing
This source code is licensed under the [GPLv2][gplv2] license. If you in any way
use this code for research that results in publications, please cite our
original article. The following [Bibtex][bibtex] entry can be used.
```
@inproceedings{Brandt2011,
  Title                    = {Wideband {MIMO} Channel Diagonalization in the Time Domain},
  Author                   = {R. Brandt and M. Bengtsson},
  Booktitle                = {Proc. IEEE Int. Symp. Indoor, Mobile Radio Commun. (PIMRC'11)},
  Year                     = {2011},
  Pages                    = {1958--1962}
}
```

## References
1. R. Brandt and M. Bengtsson, "[Wideband MIMO Channel Diagonalization in the Time Domain][pimrc2011_diva]," in _Proc. IEEE Int. Symp. Indoor, Mobile Radio Commun. (PIMRC'11)_, 2011, pp. 1958-1962.
2. J.A Foster, J.G. McWhirter, M.R. Davies, and J.A. Chambers, "[An Algorithm for Calculating the QR and Singular Value Decompositions of Polynomial Matrices][psvd]," in _IEEE Trans. Signal Process., vol.58, no.3, pp. 1263-1274, 2010
3. R. Brandt, "[Polynomial matrix decompositions: Evaluation of algorithms with an application to wideband MIMO communications][thesis_diva]," M.Sc. thesis, Uppsala University, Dec. 2010.
4. D. P. Palomar, M. A. Lagunas, A. Pascual, and A. P. Neira, "[Practical implementation of jointly designed transmit receive space-time IIR filters][palomar],” in _Proc. Int. Symp. Signal Processing and its Applications_, 2001, pp. 521-524.

[rabr5411]: http://www.kth.se/profile/rabr5411
[matben]: http://www.kth.se/profile/matben
[pimrc2011_diva]: http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-50048
[thesis_diva]: http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-134389
[psvd]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5286258
[palomar]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=950195
[matlab]: http://www.mathworks.com
[gplv2]: http://choosealicense.com/licenses/gpl-v2
[bibtex]: http://www.bibtex.org/
