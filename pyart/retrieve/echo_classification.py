"""
pyart.retrieve.echo_classification
==================================

Echo classification algorithms.

.. autosummary::
    :toctree: generated/

    echo_classification_steiner

"""

import numpy as np

from . import _echo_classification


def echo_classification_steiner(grid, refl_field, work_lev,
                                intense=40, bkg_rad=11000.0,
                                area_relation='sgp', peak_relation='sgp',
                                use_intense=True, fill_value=-9999.0):
    """
    Perform convective/statiform precipitation echo classification using
    the algorithm from Steiner et al.

    Parameters
    ----------
    grid : Grid
        Grid object on which to perform classification.
    refl_field : str
        Name of the reflectivity field in the grid.
    work_lev : float
        Elevation level, in meters, at which to perform classification.
    intense : float
        Minimum intersity, in dBZ, to automatically classify an echo as
        convective, default is 40 dBz as suggested in the Steiner article.
    bkg_rad : float
        Background radius, in meters, for calculating peakedness.  Default
        is 11,000 (11 km) as suggested in the Steiner article.
    area_relation : 'small', 'medium', 'large' or 'sgp'
        Relationship to use when determining the convective radius.  The
        'small', 'medium' and 'large' relationships correspond to the
        relationships outlined in the Steriner article.  The 'sgp'
        relationship is  an ad-hoc parameterization which has been found to
        works well for the ARM Southern Great Plains site.
    peak_relation : 'default' or 'sgp'
        Peak relationship used when calculating the peakedness metric.
        'default' uses the relationship provided in the Steiner article,
        'sgp' uses a relationship which has been found to work well for the
        ARM Southern Great Plains site.
    use_inter : bool
        True to use reflectivity intensity for classification, False to
        classify using only the peakedness metric.
    fill_val : float, optional
        Masked array fill value.

    Return
    ------
    echo_dict : dict
        Grid field dictionary containing convective/statiform echo
        classification.  Masked values indicate the gate was not
        classified, a value of 1.0 indicates convective precipitation,
        -1 stratiform precipitation.

    References
    ----------
    M. Steiner, R. A. Houze Jr., and S. E. Yuter,
    Climatological Characterization of Three-Dimensional Storm
    Structure from Operational Radar and Rain Gauge Data,
    Journal of Applied Meteorology, 1995, 34, 1978-2007.

    """
    # extract parameters needed for call to Fortran function
    x = grid.axes['x_disp']['data']
    y = grid.axes['y_disp']['data']
    z = grid.axes['z_disp']['data']
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    ze = grid.fields[refl_field]['data']
    if type(ze) is np.ma.MaskedArray:
        ze = ze.filled(fill_value)

    # perform classification
    eclass = _echo_classification.steiner(
        x, y, z, ze, dx, dy, intense, bkg_rad, work_lev, area_relation,
        peak_relation, use_intense, fill_val=fill_value)

    # create and return classification dictionary
    dic = {
        'data': np.ma.masked_equal(eclass, fill_value),
        'units': 'unitless',
        'options':
        '1: convective, -1: stratiform, masked: no classification',
        'long_name': 'Convective or stratiform precipitation classification',
    }
    return dic
