# Copyright 2017 Scott Berry, University of Zurich
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
'''Jterator module for detection of spots in images.'''
import numpy as np
import mahotas as mh
import collections
import logging
import matlab.engine

VERSION = '0.1.0'

logger = logging.getLogger(__name__)

Output = collections.namedtuple('Output', ['spots', 'figure'])


def main(image, mask, spot_size=5, rescale_quantile_min=0.01,
         rescale_quantile_max=0.99, min_of_min=np.nan, max_of_min=120,
         min_of_max=500, max_of_max=np.nan, detection_threshold=0.02,
         deblending_steps=20, plot=False, plot_clip_value=150):
    '''Detects spots as described in Battich et al., 2013. Converted to
    a jterator module from Cell Profiler by Scott Berry. Original CP
    description below.

    Parameters
    ----------
    image: numpy.ndarray[Union[numpy.uint8, numpy.uint16]]
        grayscale image in which blobs should be detected
    mask: numpy.ndarray[Union[numpy.int32, numpy.bool]]
        binary or labeled image that masks pixel regions in which spots
        should be detected
    spot_size: int, optional
        approximate spot diameter in pixels
    rescale_quantile_min: float, optional
        minimum quantile for image rescaling
    rescale_quantile_max: float, optional
        maximum quantile for image rescaling
    min_of_min: uint16, optional
        most extreme values allowed for rescaling (lowest min)
    max_of_min: uint16, optional
        most extreme values allowed for rescaling (highest min)
    min_of_max: uint16, optional
        most extreme values allowed for rescaling (lowest max)
    max_of_max: uint16, optional
        most extreme values allowed for rescaling (highest max)
    detection_threshold: uint16, optional
        threshold for pixel intensity in a spot (default: ``150``)
    deblending_steps: int, optional
        number of steps for deblending
    plot: bool, optional
        whether a plot should be generated (default: ``False``)
    plot_clip_value: bool, optional
        value to clip intensity image for plotting (default: ``150``)

    Returns
    -------
    jtmodules.identify_spots_2D.Output[Union[numpy.ndarray,numpy.ndarray,str]]

    Cell Profiler Module Help
    -------------------------

    Detects spots as destribed by Battich et al., 2013.
    ***********************************************
    Will Determine Spots in 2D Image stacks after Laplacian Of Gaussian (LoG)
    enhancing of spots. Many of the input arguments are optional. Note that
    while an external script has to be run in order to choose robust values,
    manual selection of the parameters can often yield good estimates, if
    the signal is clear enough.

    WHAT DID YOU CALL THE IMAGES YOU WANT TO PROCESS?
    Object detection should be done on this image.

    HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED PRIOR TO DEBLENDING?
    This is the name of the the spots identified after thresholding the LoG
    image.

    HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED AFTER DEBLENDING?
    Optional. Deblending can be done after spot detection to separate close
    objects. The algorithm is based upon SourceExtractor. To skip this step,
    insert / as name of the object.

    OBJECTSIZE
    This value corresponds to the approximate size of you spots. It should
    be their diameter in pixels. The LoG will use a mask of this size to
    enhance radial signal of that size. Note that in practice the specific
    value does not affect the number of spots, if spots are bright (eg.
    pixel size 5 or 6).

    INTENSITY QUANTA PER IMAGE
    Prior to spot detection the images are rescaled according to their
    intensity. Since the specific value of minimal and maximal intensities
    are frequently not robust across multiple images, intensity quantile are
    used instead. [0 1] would correspond to using the single dimmest pixel
    for minimal intensity and the single brightest pixel for maximal
    intensity. [0.01 0.90] would mean that the minimum intensity is derived
    from the pixel, which is the 1% brightest pixel of all and that the
    maximum intensity is derived from the pixel, which is the 90% brightest
    pixel.

    INTENSITY BORERS FOR INTENSITY RESCALING OF IMAGES
    Most extreme values that the image intensity minimum and image intensity
    maximum (as defined by the quanta) are allowed to have
    [LowestPossibleGreyscaleValueForImageMinimum
    HighestPossibleGreyscaleValueForImageMinimum
    LowestPossibleGreyscaleValueForImageMaximum
    HighestPossibleGreyscaleValueForImageMaximum]
    To ignore individual values, place a NaN.
    Note that these parameters very strongly depend upon the variability of
    your illumination source. When using a robust confocal microscope you can
    set the lowest and highest possible values to values,  which are very
    close (or even identical). If your light source is variable during the
    acquisition (which can be the case with Halogen lamps) you might choose
    less strict borders to detect spots of varying intensites.

    THRESHOLD OF SPOT DETECTION
    This is the threshold value for spot detection. The higher it is the more
    stringent your spot detection is. Use external script to determine a
    threshold where the spot number is robust against small variations in the
    threshold.

    HOW MANY STEPS OF DEBLENDING DO YOU WANT TO DO?
    The amount of deblending steps, which are done. The higher it is the less
    likely it is that two adjacent spots are not separated. The default of 30
    works very well (and we did not see improvement on our images with higher
    values). Note that the number of deblending steps is the main determinant
    of computational time for this module.

    WHAT IS THE MINIMAL INTENSITY OF A PIXEL WITHIN A SPOT?
    Minimal greyscale value of a pixel, which a pixel has to have in order to
    be recognized to be within a spot. Opitonal argument to make spot
    detection even more robust against very dim spots. In practice, we have
    never observed that this parameter would have any influence on the spot
    detection. However, you might include it as an additional safety measure.

    WHICH IMAGE DO YOU WANT TO USE AS A REFERENCE FOR SPOT BIAS CORRECTION?
    Here you can name a correction matrix which counteracts bias of the spot
    correction across the field of view. Note that such a correction matrix
    has to be loaded previously by a separate module, such as
    LOADSINGLEMATRIX

    Authors:
      Nico Battich
      Thomas Stoeger
      Lucas Pelkmans

    Website: http://www.imls.uzh.ch/research/pelkmans.html

    The design of this module largely follows a IdentifyPrimLoG2 by
    Baris Sumengen
    '''
    logger.debug('Parsing input variables')
    min_min = min_of_min if min_of_min > 0 else np.nan
    max_min = max_of_min if max_of_min > 0 else np.nan
    min_max = min_of_max if min_of_max > 0 else np.nan
    max_max = max_of_max if max_of_max > 0 else np.nan

    logger.debug('Starting matlab')
    mb = matlab.engine.start_matlab()
#    mb.addpath('~/matlab')
    mb.addpath('/home/tissuemaps/jtlibrary/src/matlab/cpsub/', nargout=0)

    logger.debug('Converting image to matlab format')
    image_mb = matlab.double(image.tolist())

    logger.debug('Setting options, min_of_min = %s', min_min)
    logger.debug('Setting options, max_of_min = %s', max_min)
    logger.debug('Setting options, min_of_max = %s', min_max)
    logger.debug('Setting options, max_of_max = %s', max_max)
    options = {'ObSize': float(spot_size),
               'limQuant': matlab.double([rescale_quantile_min,
                                          rescale_quantile_max]),
               'RescaleThr': matlab.double([min_min,
                                            max_min,
                                            min_max,
                                            max_max]),
               'ObjIntensityThr': matlab.uint16([]),
               'closeHoles': False,
               'ObjSizeThr': matlab.uint16([]),
               'ObjThr': detection_threshold,
               'StepNumber': deblending_steps,
               'numRatio': 0.20,
               'doLog': 0,
               'detectBias': matlab.uint16([])}

    logger.info('Detecting spots by rescaling/thresholding')
    log_filter = mb.fspecialCP3D('2D LoG', options['ObSize'])
    spots_mb = mb.ObjByFilter(image_mb,
                              log_filter,
                              options['ObjThr'],
                              options['limQuant'],
                              options['RescaleThr'],
                              options['ObjIntensityThr'],
                              True,
                              [],
                              options['detectBias'],
                              nargout=3)
    spots_pre = mb.double(mb.labelmatrix(spots_mb[1]))

    logger.info('Deblending spots using source extractor')
    spots_deblend_mb = mb.SourceExtractorDeblend(
        image_mb, spots_mb[1], spots_mb[2], options, nargout=1)

    logger.debug('Converting matlab to (contiguous) numpy array')
    # Note that there are some strange conversions between numpy and
    # matlab in terms of array byte order. Mahotas requires C arrays.
    spots = (np.array(
        spots_pre._data).reshape(
        spots_pre.size, order='F'))
    spots = np.ascontiguousarray(spots, dtype=np.int32)
    spots_deblend = (np.array(
        spots_deblend_mb._data).reshape(
        spots_deblend_mb.size, order='F'))
    spots_deblend = np.ascontiguousarray(spots_deblend, dtype=np.int32)

    logger.debug('Masking spots outside cells')
    spots[mask == 0] = 0
    mh.labeled.relabel(spots, inplace=True)
    spots_deblend[mask == 0] = 0
    mh.labeled.relabel(spots_deblend, inplace=True)

    # Exit matlab
    mb.quit()

    n_pre = spots.max()
    n_post = spots_deblend.max()
    logger.info('%d blobs detected before deblending', n_pre)
    logger.info('%d blobs detected after deblending', n_post)

    if plot:
        from jtlib import plotting

        logger.debug('dilate deblended spots')
        spots_deblend_expanded = mh.dilate(
            A=spots_deblend > 0,
            Bc=mh.disk(radius=4, dim=2))
        outlines_deblend = mh.labeled.bwperim(spots_deblend_expanded > 0)

        logger.debug('generate colorscales')
        colorscale_pre = plotting.create_colorscale(
            'Spectral', n=n_pre, permute=True, add_background=True
        )
        colorscale_post = plotting.create_colorscale(
            'Spectral', n=n_post, permute=True, add_background=True
        )

        logger.debug('create subplots')
        plots = [
            plotting.create_intensity_image_plot(
                image, 'ul',
                clip=True, clip_value=plot_clip_value
            ),
            plotting.create_mask_image_plot(
                spots, 'ur', colorscale=colorscale_pre
            ),
            plotting.create_intensity_overlay_image_plot(
                image, outlines_deblend, 'll',
                clip=True, clip_value=plot_clip_value
            ),
            plotting.create_mask_image_plot(
                spots_deblend, 'lr', colorscale=colorscale_post
            ),
        ]
        logger.info('create plot')
        figure = plotting.create_figure(
            plots,
            title='''left = before deblending ({0} spots), right = after deblending ({1} spots)'''.format(n_pre, n_post)
        )
    else:
        figure = str()

    return(Output(spots, figure))
