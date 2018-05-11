# chroma_upsampling

### Chroma upsampling for images in Y'CbCr 4:2:0 format

![Chroma upsampling](https://jarikorhonen.github.io/chromaupsample.png "Chroma upsampling")

This is a small piece of work published in ICME'15. Many practical digital images are represented in Y'CbCr 4:2:0 format, where the monochrome version (luma aka Y' component) of the image is coded in higher resolution than the color information (chroma aka Cb and Cr components). When the image is displayed, the resulting RGB image has to be reconstructed by first upsampling the low resolution Cb and Cr images to the same resolution with the Y' component. Since local variations in luma and chroma values are typically correlated, more accurate upsampling results can be obtained by exploiting the luma component in chroma upsampling process.

Matlab implementation of the chroma upsampling scheme published in ICME'18 is published here. The code should be rather self-explanatory.

If you use the implementation in your research, please cite the following publication:

J. Korhonen, "Improving image fidelity by luma-assisted chroma subsampling," *IEEE International Conference on Multimedia and Expo (ICME'15)*, Turin, Italy, June 2015. [DOI: 10.1109/ICME.2015.7177387](https://doi.org/10.1109/ICME.2015.7177387)
