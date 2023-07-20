from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import numpy as np
import copy
#from .tools import mask_good, mask_bad

import subprocess as sp 
import scipy
from scipy import stats
import inspect
import numpy as np
import healpy as hp
import os
from astropy.table import Table
import tempfile as tf
import shutil
from spt3g import core, util, mapspectra, maps, std_processing
from spt3g.mapspectra import curved_sky as cs
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from spt3g.util import healpix_tools as hpt

def average_N_spectra(spectra, N_spectra, N_ells):
    avgSpectra = np.zeros(N_ells)
    rmsSpectra = np.zeros(N_ells)
    # calcuate the average spectrum
    i = 0
    while (i < N_spectra):
        avgSpectra = avgSpectra + spectra[i, :]
        i = i + 1
    avgSpectra = avgSpectra/(1. * N_spectra)
    # calculate the rms of the spectrum
    i = 0
    while (i < N_spectra):
        rmsSpectra = rmsSpectra + (spectra[i, :] - avgSpectra)**2
        i = i + 1
    rmsSpectra = np.sqrt(rmsSpectra/(1.*N_spectra))
    #rmsSpectra = np.std(spectra, axis=0)
    return(avgSpectra, rmsSpectra)


def bin_spectrum(cls, lmin=8, lmax=None, binwidth=25, return_error=False):
    cls = np.atleast_2d(cls)
    if lmax is None:
        lmax = cls.shape[-1] - 1
    ell = np.arange(lmax + 1)
    bins = np.arange(lmin, lmax + 1, binwidth)
    ellb = stats.binned_statistic(ell, ell, statistic=np.mean, bins=bins)[0]
    clsb = np.array([stats.binned_statistic(
        ell, C, statistic=np.mean, bins=bins)[0] for C in cls]).squeeze()
    if return_error:
        clse = np.array([stats.binned_statistic(
            ell, C, statistic=np.std, bins=bins)[0] for C in cls]).squeeze()
        return ellb, clsb, clse
    return ellb, clsb

def extract_func_kwargs(func, kwargs, pop=False, others_ok=True, warn=False):
    """
    Extract arguments for a given function from a kwargs dictionary

    Arguments
    ---------
    func : function or callable
        This function's keyword arguments will be extracted.
    kwargs : dict
        Dictionary of keyword arguments from which to extract.
        NOTE: pass the kwargs dict itself, not **kwargs
    pop : bool, optional
        Whether to pop matching arguments from kwargs.
    others_ok : bool
        If False, an exception will be raised when kwargs contains keys
        that are not keyword arguments of func.
    warn : bool
        If True, a warning is issued when kwargs contains keys that are not
        keyword arguments of func.  Use with `others_ok=True`.

    Returns
    -------
    Dict of items from kwargs for which func has matching keyword arguments
    """
    spec = inspect.getargspec(func)
    func_args = set(spec.args[-len(spec.defaults):])
    ret = {}
    for k in kwargs.keys():
        if k in func_args:
            if pop:
                ret[k] = kwargs.pop(k)
            else:
                ret[k] = kwargs.get(k)
        elif not others_ok:
            msg = "Found invalid keyword argument: {}".format(k)
            raise TypeError(msg)
    if warn and kwargs:
        s = ', '.join(kwargs.keys())
        warn("Ignoring invalid keyword arguments: {}".format(s), Warning)
    return ret


def projsetup(ax=None, bm=None, projection='laea', region='spt3g', coord='C',
              lat_0=None, lon_0=None, width=None, height=None, flip='astro',
              bgcolor='0.85', graticule=True, dpar=None, dmer=None, loc_par=None,
              loc_mer=None, nsew=None, grat_opts=None, xlabel=None, ylabel=None,
              labelpad=None, xlabelpad=None, ylabelpad=None, xlabelpos=None,
              ylabelpos=None, graticule_rot=False, title='', text=None,
              text_props=None):
    """
    Create a Basemap instance and setup axes for plotting maps, including
    graticules, rotated graticules, axis labels and title text.

    Arguments
    ---------
    ax : matplotlib.Axes instance
      If None, this is assumed to be `plt.gca()`.
    bm : Basemap instance
      If None, this is created using the remaining parameters, and the
      axes instance is populated with graticules.
      If not None, the basemap instance and axes are returned without changes
      to graticules or axis labels.
    projection : string, optional
      A projection name supported by basemap. Default: equal area azimuthal
    coord : sequence of character
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    lat_0, lon_0 : float, optional
      Coordinates of map center point. Default based on region and projection.
    width, height : float, optional
      Projection plane size. Units are degree-like.
      Default based on region and projection.
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (east towards left,
      west towards right) or 'geo' (east towards right, west towards left)
    bgcolor : matplotlib color, optional
      Background colour of the image
    graticule : string, boolean, tuple, or dict, optional
      Draw a graticule in a given coordinate system. Defaults to `coord`.
      Tuple or dict behave like args to hp.graticule.
      If graticule is in different coordinates than the map, options are
      passed to `projgrat_rot` instead.
    dpar, dmer : float, optional
      Spacing of graticule parallels and meridians (degrees).
      Default based on region and projection.
    loc_par, loc_mer : 4-element list of boolean, optional
      Positions of graticule labels. [left, right, top, bottom]
      Default based on region and projection.
    grat_opts : dictionary, optional
      Dictionary of extra arguments for drawing graticules (linewidth, etc)
    nsew : bool, optional
      Whether graticule labels should use "N/S/E/W" or "+/-"
    xlabel, ylabel : string, optional
      Default based on coordinate system. Use empty string '' to turn off.
    labelpad : int, optional
      axis padding of graticule labels.  If None, use `xtick.major.pad`
      and `ytick.major.pad` from rcParams.  Override with `xoffeset`
      or `yoffset` options to `grat_opts`.
    xlabelpad, ylabelpad : int, optional
      Padding of axis labels. If None, use `axes.labelpad` from rcParams.
    xlabelpos, ylabelpos : str, optional
      Position of the axis labels: 'left', 'right', 'top' or 'bottom'
    graticule_rot : bool or dict, optional
      If True, draw rotated graticules in the alternate coordinate system
      (in G if map is in C, and vice versa).  If dict, pass arguments
      to `projgrat_rot`.
    title : str, optional
      The title of the plot. Default: ''
    text : 4 element list of str, optional
      Text to write on top of the figure in [top left, top right,
      bottom right, bottom left]. Can leave empty strings for positions
      with no text.
    text_props : dict or 4 element list of dicts, optional
      Standard text properties to apply to the text written on the figure.
      Can have different properties for the four different text locations
      if text_props is a list of dicts, or uniformly apply text_props to
      all locations if a dict.

    Returns
    -------
    ax : matplotlib.Axes instance
      The axes on which the data are to be plotted
    bm : Basemap instance
      The Basemap plotting object for the requested projection parameters
    """

    # deal with coord for region definition
    if len(coord) == 2:
        coord_region = coord[1]
    elif len(coord) == 1:
        coord_region = coord
    else:
        raise ValueError("Invalid 'coord' argument")
    if coord_region not in ['C', 'G']:
        raise ValueError("Invalid 'coord' argument")

    # select axes labels
    xlabeldict = {'G': 'Galactic Longitude', 'C': 'Right Ascension'}
    ylabeldict = {'G': 'Galactic Latitude', 'C': 'Declination'}

    # graticule options
    graticule = coord_region if graticule is True else graticule
    if isinstance(graticule, tuple):
        dpar = graticule[0]
        dmer = graticule[1]
        graticule = graticule[2]
    elif isinstance(graticule, dict):
        graticule = graticule.copy()
        dpar = graticule.pop("dpar", None)
        dmer = graticule.pop("dmer", None)
        grat_opts = graticule if grat_opts is None else grat_opts
        graticule = graticule.pop("coord", coord_region)
    if grat_opts is None:
        grat_opts = {}
    else:
        grat_opts = grat_opts.copy()

    if projection == 'moll':
        region = None
        width = None
        height = None
        dpar = 30 if dpar is None else dpar
        dmer = 30 if dmer is None else dmer
        loc_par = [0, 0, 0, 0] if loc_par is None else loc_par
        loc_mer = [0, 0, 0, 0] if loc_mer is None else loc_mer
        lon_0 = 0 if lon_0 is None else lon_0
    elif projection == 'laea':
        loc_par = [1, 0, 0, 1] if loc_par is None else loc_par
        loc_mer = [0, 0, 1, 0] if loc_mer is None else loc_mer
        xlabel = xlabeldict[coord_region] if xlabel is None else xlabel
        ylabel = ylabeldict[coord_region] if ylabel is None else ylabel
        xlabelpos = 'top' if xlabelpos is None else xlabelpos
        ylabelpos = 'left' if ylabelpos is None else ylabelpos
        xlabelpad = 25 if xlabelpad is None else xlabelpad
        ylabelpad = 30 if ylabelpad is None else ylabelpad
        dpar = 15 if dpar is None else dpar
        dmer = 20 if dmer is None else dmer
        if region in ['spt3g','spt3g_summer']:
            if coord_region == 'C':
                lat_0 = -60 if lat_0 is None else lat_0
                lon_0 = 0 if lon_0 is None else lon_0
                width = 80 if width is None else width
                height = 40 if height is None else height
            elif coord_region == 'G':
                lat_0 = -75 if lat_0 is None else lat_0
                lon_0 = -60 if lon_0 is None else lon_0
                width = 120 if width is None else width
                height = 40 if height is None else height
        elif region == 'cmb_mask':
            if coord_region == 'C':
                lat_0 = -37 if lat_0 is None else lat_0
                lon_0 = 52 if lon_0 is None else lon_0
                width = 72 if width is None else width
                height = 48 if height is None else height
            elif coord_region == 'G':
                lat_0 = -56 if lat_0 is None else lat_0
                lon_0 = -118 if lon_0 is None else lon_0
                width = 50 if width is None else width
                height = 72 if height is None else height
        elif region is not None:
            raise ValueError("Invalid region: {} for {}".format(
                region, projection))
    else:
        warn("Guessing defaults for projection {}. May need to set"
             " lat_0, lon_0, width, height, etc")
        lat_0 = -50 if lat_0 is None else lat_0
        lon_0 = 65 if lon_0 is None else lon_0
        width = 135 if width is None else width
        height = 90 if height is None else height
        dpar = 30 if dpar is None else dpar
        dmer = 30 if dmer is None else dmer

    # rotated graticule options
    if graticule and graticule != coord_region:
        graticule = False
        graticule_rot = grat_opts
        if dpar is not None:
            graticule_rot.setdefault('dpar', dpar)
        if dmer is not None:
            graticule_rot.setdefault('dmer', dmer)
    if isinstance(graticule_rot, dict):
        grat_opts_rot = graticule_rot
        graticule_rot = True
    else:
        grat_opts_rot = {}
    if projection == 'moll' and graticule_rot:
        # buggy
        NotImplementedError(
            "Different grat coords not implemented for mollweide projection")

    if graticule:
        style = "+/-" if not nsew else ""
        grat_opts.setdefault('dpar', dpar)
        grat_opts.setdefault('dmer', dmer)
        grat_opts.setdefault('loc_par', loc_par)
        grat_opts.setdefault('loc_mer', loc_mer)
        grat_opts.setdefault('flip', flip)
        grat_opts.setdefault('fontsize', 'medium')
        grat_opts.setdefault('labelstyle', style)
        grat_opts.setdefault('labelpad', labelpad)
        grat_opts.setdefault('label_par', ylabel)
        grat_opts.setdefault('label_mer', xlabel)
        grat_opts.setdefault('labelpos_par', ylabelpos)
        grat_opts.setdefault('labelpos_mer', xlabelpos)
        grat_opts.setdefault('labelpad_par', ylabelpad)
        grat_opts.setdefault('labelpad_mer', xlabelpad)

    if graticule_rot:
        grat_opts_rot.setdefault(
            'coord', ['G' if coord_region == 'C' else 'C', coord_region])

    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    if ax is None:
        ax = plt.gca()

    # create projection and set up axes
    if bm is None:
        bm = Basemap(projection=projection, rsphere=180/np.pi, lat_0=lat_0,
                     lon_0=lon_0, width=width, height=height, resolution=None)

        # setup axes
        bm.set_axes_limits(ax=ax)
        if flip == "astro":
            xl = ax.get_xlim()
            if xl[0] < xl[1]:
                # fix and invert axis limits
                ax.set_xlim(xl[1], xl[0])
        if bgcolor:
            ax.set_facecolor(bgcolor)

        # draw graticules
        if graticule:
            projgrat(bm, ax=ax, **grat_opts)
        if graticule_rot:
            projgrat_rot(bm, ax=ax, **grat_opts_rot)

        # store view configuration
        bm.region = region
        bm.coord = coord_region

    # figure text
    if title:
        plt.suptitle(title)
    if text_props is None:
        text_props = {}
    if text is not None:
        if isinstance(text, str):
            text = [text, "", "", ""]
        textlocs = {0: ['top', 'left', 0.02, .98],
                    1: ['top', 'right', 0.98, 0.98],
                    2: ['bottom', 'right', 0.98, 0.02],
                    3: ['bottom', 'left', 0.02, 0.02]}
        if isinstance(text_props, dict):
            text_props = text_props.copy()
            text_props = np.tile(text_props, 4)
        else:
            for i in xrange(len(text_props)):
                text_props[i] = text_props[i].copy()
        for l in range(4):
            text_props[l].setdefault('zorder', 5)
            ax.text(textlocs[l][2], textlocs[l][3], text[l],
                    verticalalignment=textlocs[l][0],
                    horizontalalignment=textlocs[l][1],
                    transform=ax.transAxes, **text_props[l])

    # return
    return ax, bm


def projview(m, vmin=None, vmax=None, pmin=0.5, pmax=99.5, coord='C',
             resol=None, diverging=True, cmap=None, log=False, cbar=True, unit=None,
             cbar_extend='neither', cbar_ticks=None, cbar_minorticks=False,
             nest=False, return_projected_map=False, interp=True, pad=0.0, **kwargs):
    """
    Plot an healpix map (given as an array) in arbitrary projection.
    Uses mpl_toolkits.basemap to project map and graticules.
    Options are similar to healpy.cartview, with a few adjustments for
    SPT3G-specific plotting (e.g. the 'region' option) and projections.

    Parameters
    ----------
    map : float, array-like or None
      An array containing the map,
      supports masked maps, see the `ma` function.
      If None, will display a blank map, useful for overplotting.
    pmin, pmax : float, optional
      Percentile minimum and maximum range values. Can also provide as a
      string with "%" to (v)min or (v)max.
    min, max, vmin, vmax : float, optional
      The minimum and maximum range values. Override pmin, pmax.
      v-less forms for compatibility with healpy.
    coord : sequence of character
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
      If the coordinate system of the map differs from that of the
      configured `Basemap`, the map coordinates will be rotated.
    resol : float, optional
      Resolution of the image in pixels per projected "degree". Image will
      have dimenstion (width * resol) by (height * resol).
    diverging : bool, optional
      Whether to set symmetric colour limits for diverging colour scale.
      Default True.
    cmap : string or colormap, optional
      The colormap to use. Default chosen based on region
    log : bool, optional
    norm : `matplotlib.colors.Normalize` instance
      If `log` is True, use a logarithmic colour scale, otherwise linear.
      For other general normalizations use the norm kwarg, which will be
      passed to imshow (and override this option)
    cbar : bool or string, optional
      Display the colorbar. String specifies position ('left', 'right',
      'bottom', 'top'). Will interpret 'vertical' as 'right' or 'horizontal'
      as 'bottom'. Default: 'right' when set to True.
    unit : str, optional
      A text describing the unit of the data.
    cbar_extend : string, optional
      Set to "min", "max", or "both" to show extended colour for out-of-range
    cbar_ticks : list, optional
      Colorbar tick positions.
    cbar_minorticks : bool, optional
      If True, turn on minor ticks on colorbar axis.
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    return_projected_map : bool, optional
      if True returns the projected map in a 2d numpy array
    interp : bool, optional
      If True, interpolate values using bilinear interpolation.  If False,
      use the nearest neighbor value.

    Keyword Arguments:
    ------------------
    Remaining arguments are passed to `projsetup` to create the Basemap
    instance and populate the figure axes with labels and graticules.
    Notably:

    bm : Basemap instance
      A pre-configured Basemap instance
    projection : string, optional
      A projection name supported by basemap. Default: equal area azimuthal
    region : string, optional
      Determine projection parameters automatically based on the desired
      region to plot. Options are 'spider' (SPIDER CMB region aka 'cmb')
      and 'rcw38'. Can specify 'both' for CMB + RCW38.
    coord : sequence of character
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.

    Any further remaining arguments are passed to imshow. Notably:

    ax : matplotlib Axes instance
        The axes in which to plot. For use with, eg, pyplot.subplots

    Returns
    -------
    mxy : array_like, optional
      Returned if return_projected_map is True
    im : matplotlib.image.AxesImage instance
      The image object
    bm : basemap.Basemap
      The basemap object. For extra plotting functions in map coords
    cbar : matplotlib.colorbar.Colorbar instance, optional
      The colorbar object, if cbar is set
    """

    import matplotlib.pyplot as plt

    setup_opts = extract_func_kwargs(
        projsetup, kwargs, pop=True, others_ok=True)
    ax, bm = projsetup(coord=coord, **setup_opts)
    projection = bm.projection
    region = bm.region
    if list(coord)[-1] != bm.coord:
        coord = [list(coord)[0], bm.coord]

    if cmap is None:
        cmap = plt.cm.RdBu_r if diverging else "inferno"

    # colorbar options
    if cbar is True or cbar in ['v', 'vertical', 'r']:
        cbar = 'right'
    elif cbar in ['h', 'horizontal', 'b','bottom']:
        cbar = 'bottom'
        pad = 0.4
    cbar_label = unit if unit is not None else ""

    # color scale
    vmin = kwargs.pop('min', None) if vmin is None else vmin
    vmax = kwargs.pop('max', None) if vmax is None else vmax
    if isinstance(vmin, str) and '%' in vmin:
        pmin = float(vmin.rstrip('%'))
    if isinstance(vmax, str) and '%' in vmax:
        pmax = float(vmax.rstrip('%'))
    if pmin is not None or pmax is not None:
        mask = mask_good(m, badinf=True)
        if vmin is None:
            vmin = np.percentile(m[mask], pmin)
        if vmax is None:
            vmax = np.percentile(m[mask], pmax)

    if diverging:
        val = np.max([np.abs(vmin), np.abs(vmax)])
        vmin = -val
        vmax = val

    if kwargs.get("norm", None) == "log":
        kwargs.pop("norm")
        log = True
    if log:
        from matplotlib.colors import LogNorm
        kwargs.setdefault("norm", LogNorm(vmin=vmin, vmax=vmax))

    # convert UNSEEN values to NAN (temorarily)
    m = m.copy()
    unseen_mask = (m == hp.UNSEEN)
    m[unseen_mask] = np.nan
    if resol is None:
        resol = 5 if projection == 'moll' else 15
    nx = int((bm.xmax - bm.xmin) * resol)
    ny = int((bm.ymax - bm.ymin) * resol)
    phi, theta = bm.makegrid(nx, ny)
    # convert lat/long from makegrid to theta/phi in radians
    theta = np.deg2rad(90 - theta)
    phi = np.deg2rad(phi)

    # rotate pixel angles if converting coordinates
    if len(coord) == 2:
        # note that I need rotation for opposite direction from coord
        rmat = hp.rotator.get_coordconv_matrix((coord[1], coord[0]))[0]
        theta, phi = hp.rotator.rotateDirection(
            rmat, theta.ravel(), phi.ravel())
        theta = theta.reshape((ny, nx))
        phi = phi.reshape((ny, nx))

    # get values at each projected pixel. Either bilinear or nearest neighbour
    if interp:
        dat = hp.get_interp_val(m, theta, phi, nest=nest)
    else:
        dat = m[hp.ang2pix(hp.npix2nside(len(m)), theta, phi, nest=nest)]

    im = bm.imshow(dat, vmin=vmin, vmax=vmax, origin="lower", cmap=cmap,
                   ax=ax, **kwargs)
    if projection == 'moll':
        limb = bm.drawmapboundary(fill_color='white', ax=ax)
        im.set_clip_path(limb)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
    
    if cbar:
        cb = bm.colorbar(im, location=cbar, label=cbar_label,
                         extend=cbar_extend, pad=pad, ax=ax, format='%.0e')
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
        if cbar_minorticks:
            cb.ax.minorticks_on()
        plt.sca(ax)     # reset current axis to main image
    else:
        cb = None

    # restore UNSEEN pixels
    m[unseen_mask] = hp.UNSEEN

    plt.draw()

    # return stuff
    return (dat,) * return_projected_map + (im, bm) + (cb,) * bool(cbar)


def projgrat(bmap, dpar=15, dmer=20, loc_par=[1, 0, 0, 1], loc_mer=[0, 0, 1, 0],
             label_par=None, label_mer=None, labelpad=None, flip='astro',
             zorder=5, ax=None, **kwargs):
    '''Draw lat-lon graticules.

    In Basemap, the label pad is computed in projection units. Now you can use
    the keyword argument 'labelpad' to control this separation in points. If
    not specified then this value is taken from rcParams.

    Arguments:

    bmap : Basemap object
    dpar, dmer : int
        Difference in degrees from one longitude or latitude to the next.
    loc_par, loc_mer : 4-tuple of bools
        list of 4 values (default [0,0,0,0]) that control
        whether parallels are labelled where they intersect
        the left, right, top or bottom of the plot. For
        example labels=[1,0,0,1] will cause parallels
        to be labelled where they intersect the left and
        and bottom of the plot, but not the right and top
    label_par, label_mer : string, optional
        Parallel and meridian axis labels.
    labelpad : int, optional
        axis padding of graticule labels.  If None, use `xtick.major.pad`
        and `ytick.major.pad` from rcParams.
    labelpad_par, labelpad_mer : int, optional
        axis padding of axis labels.  If None, use `axes.labelpad` from
        rcParams.
    labelpos_par, labelpos_mer : string, optional
        location of axis labels. 'left', 'right', 'top' or 'bottom'.
        if not supplied, these are automatically determined based on `loc_par`
        and `loc_mer`.
    flip : {'astro', 'geo'}, optional
        Defines the convention of projection : 'astro' (east towards left,
        west towards right) or 'geo' (east towards right, west towards left)
    ax : matplotlib Axes instance
        Axes in which to draw graticules. Current axes by default.
    **kwargs -- Other arguments to drawparallels, drawmeridians and plt.text.
    '''
    # Processes arguments and rcParams for default values
    import matplotlib.pyplot as plt
    # kwargs.setdefault('color', plt.rcParams['grid.color'])
    # kwargs.setdefault('linewidth', plt.rcParams['grid.linewidth'])

    if ax is None:
        ax = plt.gca()

    if labelpad is not None:
        padx = pady = labelpad
    else:
        pady = plt.rcParams['xtick.major.pad']
        padx = plt.rcParams['ytick.major.pad']
    labelsize = kwargs.pop('fontsize', None)
    if labelsize is not None:
        xfontsize = yfontsize = labelsize
    else:
        xfontsize = plt.rcParams['xtick.labelsize']
        yfontsize = plt.rcParams['ytick.labelsize']

    labelpad_par = kwargs.pop('labelpad_par', None)
    if labelpad_par is None:
        labelpad_par = plt.rcParams['axes.labelpad']
    labelpad_mer = kwargs.pop('labelpad_mer', None)
    if labelpad_mer is None:
        labelpad_mer = plt.rcParams['axes.labelpad']
    labelpos_par = kwargs.pop('labelpos_par', None)
    labelpos_mer = kwargs.pop('labelpos_mer', None)

    # Vectors of coordinates
    parallels = np.arange(0., 90., dpar)
    parallels = np.sort(np.concatenate((parallels, -parallels[1:])))
    meridians = np.arange(0., 180., dmer)
    meridians = np.sort(np.concatenate((meridians, -meridians[1:])))
    kwargs.setdefault(
        'latmax', parallels[-1] if parallels[-1] <= 85 else parallels[-2])

    if loc_par is None:
        loc_par = [0] * 4
    if loc_mer is None:
        loc_mer = [0] * 4
    if labelpos_par is None:
        labelpos_par = ('left' if loc_par[0] else 'right' if loc_par[1]
                        else 'bottom' if loc_par[3] else 'top' if loc_par[2]
                        else 'left')
    if labelpos_mer is None:
        labelpos_mer = ('bottom' if loc_mer[3] else 'top' if loc_mer[2]
                        else 'left' if loc_mer[0] else 'right' if loc_mer[1]
                        else 'bottom')

    if flip == "astro":
        # also need to flip graticule labels
        loc_par[0], loc_par[1] = loc_par[1], loc_par[0]
        loc_mer[0], loc_mer[1] = loc_mer[1], loc_mer[0]

    # If not specified then compute the label offset by 'labelpad'
    xos = kwargs.pop('xoffset', None)
    yos = kwargs.pop('yoffset', None)
    if xos is None and yos is None:
        # Page size in inches and axes limits
        fig_w, fig_h = plt.gcf().get_size_inches()
        (x1, y1), (x2, y2) = plt.gca().get_position().get_points()
        # Width and height of axes in points
        w = (x2 - x1) * fig_w * 72
        h = (y2 - y1) * fig_h * 72
        # If the aspect relation is fixed then compute the real values
        if bmap.fix_aspect:
            aspect = bmap.aspect * w / h
            if aspect > 1:
                w = h / bmap.aspect
            elif aspect < 1:
                h = w * bmap.aspect
        # Offset in projection units (meters or degrees)
        xos = padx * (bmap.urcrnrx - bmap.llcrnrx) / w
        yos = pady * (bmap.urcrnry - bmap.llcrnry) / h

    # Draws the grid
    retpar = bmap.drawparallels(parallels, labels=loc_par, fontsize=yfontsize,
                                xoffset=xos, yoffset=yos, ax=ax, **kwargs)
    if flip == 'astro':
        for ll, tt in retpar.values():
            for t in tt:
                if t.get_ha() == 'left':
                    t.set_ha('right')
                elif t.get_ha() == 'right':
                    t.set_ha('left')
    bmap.parlabels = retpar

    retmer = bmap.drawmeridians(meridians, labels=loc_mer, fontsize=xfontsize,
                                xoffset=xos, yoffset=yos, ax=ax, **kwargs)
    if flip == 'astro':
        for ll, tt in retmer.values():
            for t in tt:
                if t.get_ha() == 'left':
                    t.set_ha('right')
                elif t.get_ha() == 'right':
                    t.set_ha('left')
    bmap.merlabels = retmer

    def draw_label(label, labelpos, labelpad):
        if labelpos in ['left', 'right']:
            ax.set_ylabel(label, labelpad=labelpad)
            ax.yaxis.set_label_position(labelpos)
        else:
            ax.set_xlabel(label, labelpad=labelpad)
            ax.xaxis.set_label_position(labelpos)

    if label_par:
        draw_label(label_par, labelpos_par, labelpad_par)
    if label_mer:
        draw_label(label_mer, labelpos_mer, labelpad_mer)


def mask_good(m, badval=hp.UNSEEN, rtol=1.e-5, atol=1.e-8,
              badnan=True, badinf=True):
    """Returns a bool array with ``False`` where m is close to badval,
    NaN or inf.

    Parameters
    ----------
    m : a map (may be a sequence of maps)
    badval : float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol : float, optional
        The relative tolerance
    atol : float, optional
        The absolute tolerance
    badnan : bool, optional
        If True, also mask NaN values
    badinf : bool, optional
        If True, also mask inf values

    Returns
    -------
    a bool array with the same shape as the input map, ``False`` where input map is
    close to badval, NaN or inf, and ``True`` elsewhere.

    See Also
    --------
    mask_bad, ma

    Examples
    --------
    >>> import healpy as hp
    >>> m = np.arange(12.)
    >>> m[3] = hp.UNSEEN
    >>> m[4] = np.nan
    >>> mask_good(m)
    array([ True,  True,  True, False, False,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)
    """
    m = np.asarray(m)
    mask = np.ones_like(m, dtype=bool)
    if badnan:
        mask &= ~np.isnan(m)
    if badinf:
        mask &= np.isfinite(m)
    mask[mask] = hp.mask_good(m[mask], badval=badval, rtol=rtol, atol=atol)
    return mask

def fit_foreground_template(m, template, p0, mask=None, joint_fit=True,
    smooth=False, fwhm=.5, in_place=False):
    '''
    Fits the template to the input map by minimizing the variance in the
    cleaned map

    Args
    ----------
    m : float array
        The input map.
    template : list of float_array
        The foreground templates. It should be leakage subtracted.
    p0 : list of float
        The initial guess for the fit.

    Optional Args
    -------------
    mask : bool array
        The region to mask. Default None.
    joint_fit : bool
        Fit I and P at the same time. Default True.
    smooth : bool
        Presmooth the map before fitting. Smooths both m and template. Default
        is False
    fwhm : float
        If smooth is True, then smooth with a beam with full-width half max
        fwhm in arcmin. Default is 1.2, this is the one of SPT3G 150GHz.
    in_place : bool
        If True, does not create a copy of the data to perform operations. This
        will alter the input maps.


    Returns
    -------
    fit_amp : float
        Returns the fit parameter. If joint_fit is False, returns two fit
        values, the first for I and the second for P. If the fit fails,
        returns np.nan
    '''
    m = np.array([np.asarray(m['T']), np.asarray(m['Q']), np.asarray(m['U'])])
    template = np.array([np.asarray(template['T']), np.asarray(template['Q']), np.asarray(template['U'])])

    if mask is not None:
        if not isinstance(mask[1], np.bool_):
            mask = ~mask.astype(bool)
        else:
            mask = ~mask
        m[..., mask] = np.nan

    if smooth:
       m = hpt.smoothing(m, fwhm=np.deg2rad(fwhm))

    # Apply mask
    if mask is not None:
        template[..., mask] = np.nan
    
    if smooth:
        template[:] = hpt.smoothing(template, fwhm=np.deg2rad(fwhm))

    # Define the residual
    def residual(alpha, m, templates, mode):
        dm = np.copy(m)
    
        dm -= alpha * templates
    
        map_variance = np.ones(len(m))    
        for i, d in enumerate(dm):
            if mask is not None:
                map_variance[i] = np.nansum((d[~mask] - np.nanmean(d[~mask]))**2)
            else:
                map_variance[i] = np.nansum((d - np.nanmean(d))**2)
        if mode == 't':
            return map_variance[0]
        elif mode == 'p':
            return map_variance[1] + map_variance[2]
        else:
            return map_variance.sum()  

    # Solve for the fit parameter
    from scipy.optimize import minimize
    if joint_fit:
        f = minimize(residual, p0, args=(m, template, 1), method='Nelder-Mead')
        fit_amp = np.nan
        if f.success:
            fit_amp = f.x 
    else:
        f_I = minimize(residual, p0, args=(m, template, 't'), method='Nelder-Mead')
        f_P = minimize(residual, p0, args=(m, template,'p' ), method='Nelder-Mead')
        fit_amp = np.zeros(2) * np.nan
        if f_I.success:
            fit_amp[0] = f_I.x
        if f_P.success:
            fit_amp[1] = f_P.x
    return fit_amp

def pts_mask_hfi(mmap, freq, input_coord, nside=2048, n=None):
    """
    For HFI maps only. This function masks the n brighest (psfflux) points in mmap
    with radius given by the planck beam.

    Arguments
    =========
    mmap : numpy array
        The input map to get point source masked.
    freq : int
        The pfreq channel to pick point sources from.
    input_coord : string
        Input coord. If not 'C', then the maps will be rotated to 'C', point sources are
        masked, then rotate back to 'G' in the output
    nside : int, optional
        The nside of the input map.
    n : int, optional
        The number of brightest point source to be masked.

    Returns
    =======
    The map with point source masked.
    """
    
    if input_coord == 'G':
        try:
            mmap = rotate_map(mmap, coord=('G', 'C'))
        except:
            mmap = rotate_map(GtoIQU(mmap), coord=('G', 'C'))
    path_to_points=os.path.join('/sptgrid/user/rgualtie/PlanckDR3Maps','planck_point_sources/')
    pointspath=os.listdir(path_to_points)
    points={}

    for j in pointspath:
        if '.tbl' in j:
            key=int(j.strip('IpacTableFromSource_.t'))
            points[key]=Table.read(os.path.join(path_to_points,j), format='ipac')

    if not isinstance(freq, int):
        freq = int(freq)
    if not isinstance(n, int):
        n = int(n)

    for i in points:
        points[i]=points[i][(points[i]['ra']>10) & (points[i]['ra']<100) & \
            (points[i]['dec']<-16) & (points[i]['dec']>-58)]
        points[i].sort('psfflux')
        points[i].reverse()
    r = 0.03*3 #100 arcmin radius

    try:
        cutmap=np.array(mmap)
    except:
        cutmap=np.array(G3toIQU(mmap))
    for j,i in enumerate(points[freq][:n]):
        ra=i['ra']
        dec=i['dec']
        theta=-dec*np.pi/180+np.pi/2
        phi=ra*np.pi/180
        vec=hp.ang2vec(theta,phi)
        r=0.2*np.pi/180
        disc=hp.query_disc(nside,vec,r)
        for k in cutmap:
            k[disc]=hp.UNSEEN
            while (k[disc]==hp.UNSEEN).any():
                for ii in disc:
                    if k[ii]==hp.UNSEEN:
                        neighbs=hp.pixelfunc.get_all_neighbours(nside,ii)
                        k[ii]=np.mean(k[neighbs][k[neighbs]!=hp.UNSEEN])
    if input_coord == 'G':
        cutmap = tools.rotate_map(cutmap, coord=('C', 'G'))

    return cutmap

def ell2deg(ell):
    return 180./60.*np.sqrt(2./np.pi/ell)

def plot_spectra(mm1=None, mm2=None, lowell=True, title='', label=''):
    '''Input a map or two maps and plot the spectra'''
    ellb, dlsb, errb = cs.spectrum_spice(map1=mm1, map2=mm2, lmin=8, lmax=1500, bin_width=25, mask=mask, lfac=True, return_error=True)
    fig, ax = plt.subplots(2,3, sharex=True)
    titles = ['TT','EE','BB','TE','EB','TB']
    plt.suptitle(title)
    for i in range(6):
        ax[i//3,i%3].set_title(titles[i])
        ax[i//3,i%3].plot(th_specs[i], 'k-', label='$\Lambda$-CDM')
        ax[i//3,i%3].errorbar(ellb, dlsb[i], yerr=errb[i], fmt='o', label=label)
        if i>2:
            ax[i//3,i%3].set_xlabel('Multipole $\ell$')
        if i==0 or i==3:
            ax[i//3,i%3].set_ylabel('$D_\ell$ $[\mu K^2]$')
        if lowell:
            ax[i//3,i%3].set_xlim(0,500)
        else:
            ax[i//3,i%3].set_xlim(0,1500)
        ax[i//3,i%3].legend(); ax[i//3,i%3].grid()
    return True


"""
Tools for FITS world coordinte system 2D projections. Mainly for plotting.

Requires astropy.

For more info on WCS, see: https://arxiv.org/pdf/astro-ph/0207413.pdf
"""

__all__ = ["get_wcs_for_region", "healpix_to_wcs", "wcs_subplot",
        "wcs_plot_world", "wcsview"]

def get_wcs_for_region(region="cmb", coord="C", projection="ZEA", **kwargs):
    """
    Get the pre-defined WCS object for a named region

    Arguments
    =========
    region : str
        The name of the region. Default: "cmb"
    coord : str
        Coordinate system: "C" for equatorial (the default),
        or "G" for galactic, "E" for ecliptic
    projection : str
        The WCS code for the projection type.
        Defaults to "ZEA": zenithal (azimuthal) equal area.

    Keyword arguments are used to update WCS parameters from the defaults.
    For example: crval=[42., -42.]

    Returns
    =======
    an astropy.wcs.WCS object for the region
    """
    from astropy import wcs

    wcs_params = _get_wcs_params_for_region(region=region, coord=coord,
            projection=projection)
    wcs_params.update(kwargs)

    # naxis needs to be set at construction time, in header
    naxis = wcs_params.pop("naxis")
    header = {"NAXIS{}".format(i+1): n for (i, n) in enumerate(naxis)}

    w = wcs.WCS(header=header, fix=False)
    for attr, val in wcs_params.items():
        setattr(w.wcs, attr, val)

    # store region name as a non-standard attribute
    w.region = region
    return w

def _get_wcs_params_for_region(region="cmb", coord="C", projection="ZEA"):
    """
    Helper function that is a big lookup of WCS parameters for region.

    See get_wcs_for_region for more info and usage.
    """
    # handle ctype first
    if coord == "C":
        cnames = ["RA--", "DEC-"]
    elif coord in ["G", "E"]:
        cnames = [n.format(coord) for n in ["{}LON", "{}LAT"]]
    else:
        raise ValueError("Unknown coordinate system: {}".format(coord))
    ctype = ["{}-{}".format(c, projection) for c in cnames]

    # lookup based on (region, coord, projection)
    # Explanation for some parameters (requirements vary between projections):
    #     naxis: number of pixels along each axis
    #     crval: world/sky coordinates of reference pixel
    #     crpix: image/pixel coordinates of reference pixel
    #     cdelt: pixel size at reference pixel
    lookup = {
        # ZEA = Zenithal (azimuthal) equal area
        ("cmb", "C", "ZEA"): dict(naxis=[800, 400],
                #crval=[50., -40.], crpix=[435, 315], cdelt=[-0.1, 0.1]),
                crval=[0., -50.], crpix=[400, 300], cdelt=[-0.1, 0.1]),
        ("latlon", "C", "ZEA"): dict(naxis=[620, 460],
                crval=[49., -36.], crpix=[310, 235], cdelt=[-0.1, 0.1]),
        ("rcw38", "C", "ZEA"): dict(naxis=[550, 420],
                crval=[131., -49.], crpix=[275, 210], cdelt=[-0.1, 0.1]),
        ("both", "C", "ZEA"): dict(naxis=[1350, 900],
                crval=[65., -50.], crpix=[675, 450], cdelt=[-0.1, 0.1]),
        ("cmb", "G", "ZEA"): dict(naxis=[655, 970],
                crval=[-127., -63.], crpix=[380, 485], cdelt=[-0.1, 0.1]),
        ("latlon", "G", "ZEA"): dict(naxis=[460, 670],
                crval=[-129., -56.], crpix=[270, 360], cdelt=[-0.1, 0.1]),
        ("rcw38", "G", "ZEA"): dict(naxis=[550, 550],
                crval=[-93., -2.], crpix=[275, 275], cdelt=[-0.1, 0.1]),
    }
    # defaults, independent of region, just (coord, projection)
    defaults = {
        # ZEA = Zenithal (azimuthal) equal area
        ("C", "ZEA"): dict(naxis=[970, 650],
                crval=[50., -40.], crpix=[435, 315], cdelt=[-0.1, 0.1]),
        ("G", "ZEA"): dict(naxis=[655, 970],
                crval=[-127., -63.], crpix=[380, 485], cdelt=[-0.1, 0.1]),
        # MOL = Mollweide
        ("C", "MOL"): dict(naxis=[650, 326],
                crval=[0., 0.], crpix=[325, 163], cdelt=[-0.5, 0.5]),
        ("G", "MOL"): dict(naxis=[650, 326],
                crval=[0., 0.], crpix=[325, 163], cdelt=[-0.5, 0.5]),
    }

    key = (region, coord, projection)
    try:
        params = lookup[key]
    except KeyError:
        try:
            params = defaults[(coord, projection)]
        except KeyError:
            raise ValueError("Unsupported (region, coord, projection): {}".format(key))
    params["ctype"] = ctype
    return params

def _get_grid_spacing_for_region(region="cmb", coord="C", projection="ZEA"):
    """
    Helper function that is a big lookup of grid spacing.

    Returns
    =======
    (x_space, y_space) in degrees
    """
    # region specific values
    lookup = {
        ("latlon", "C", "ZEA"): (20., 12.),
    }
    # defaults, independent of region, just (coord, projection)
    defaults = {
        ("C", "ZEA"): (40., 15.),
        ("G", "ZEA"): (45., 15.),
        ("C", "MOL"): (60., 30.),
        ("G", "MOL"): (60., 30.),
    }
    try:
        return lookup[(region, coord, projection)]
    except KeyError:
        try:
            return defaults[(coord, projection)]
        except KeyError:
            return None

def healpix_to_wcs(hmap, w):
    """
    Interpolate a healpix map to a WCS grid

    Arguments
    =========
    hmap : array
        Map in healpix pixelization
    w : astropy.wcs.WCS
        WCS pixelization information. Eg from get_wcs_for_region
    """
    import healpy as hp
    npixx, npixy = w._naxis
    x, y = np.arange(npixx), np.arange(npixy)
    x, y = np.meshgrid(x, y)
    phi, theta = w.all_pix2world(x, y, 0)
    phi = np.radians(phi)
    theta = np.pi/2 - np.radians(theta)
    # some projections have nonsense values in some places
    mask = np.logical_and(np.isfinite(phi), np.isfinite(theta))
    img = np.full_like(phi, np.nan)
    img[mask] = hp.get_interp_val(hmap, theta[mask], phi[mask])
    return img

def wcs_subplot(w, grid_space=None, grid_opts=None, subplot=(), **kwargs):
    """
    Wrapper around pyplot.subplot that sets WCS information

    Arguments
    =========
    w : astropy.wcs.WCS
        WCS pixelization information. Eg from get_wcs_for_region
    grid_space : (float, float), optional
        The spacing of (x, y) grid lines, in degrees.
        Will guess sensible defaults if not specified.
    grid_opts : dict, optional
        Dictionary of style options for the coordinate grid
    subplot : tuple, optional
        Positional argument typle to pass to pyplot.subplot.
        Eg (nrow, ncol, index)

    Remaining positional and keyword arguments are passed pyplot.subplot

    Returns
    =======
    the subplot axes
    """
    import matplotlib.pyplot as plt
    import astropy.units as u

    # check for projections that should use an elliptical frame
    projection = w.wcs.ctype[0][-3:]
    if projection in ["MOL"]:
        from astropy.visualization.wcsaxes.frame import EllipticalFrame
        kwargs["frame_class"] = EllipticalFrame

    kwargs["projection"] = w
    ax = plt.subplot(*subplot, **kwargs)

    # change some of the default axis grid configuration
    ax.coords[0].set_major_formatter("dd")
    ax.coords[1].set_major_formatter("dd")
    if w.wcs.ctype[0].startswith("RA"):
        coord = "C"
        ax.coords[0].set_axislabel("Right Ascension")
        ax.coords[1].set_axislabel("Declination")
    elif w.wcs.ctype[0].startswith("GLON"):
        coord = "G"
        ax.coords[0].set_axislabel("Galactic Longitude")
        ax.coords[1].set_axislabel("Galactic Latitude")
    elif w.wcs.ctype[0].startswith("ELON"):
        coord= "E"
        ax.coords[0].set_axislabel("Ecliptic Longitude")
        ax.coords[1].set_axislabel("Ecliptic Latitude")
    else:
        coord = None

    # check for region-based grid spacing defaults
    if grid_space is None:
        try:
            region = getattr(w, "region")
        except AttributeError:
            region = None
        grid_space = _get_grid_spacing_for_region(region, coord, projection)
    if grid_space is not None:
        ax.coords[0].set_ticks(spacing=grid_space[0]*u.deg)
        ax.coords[1].set_ticks(spacing=grid_space[1]*u.deg)

    # set default grid line options
    if grid_opts is None:
        grid_opts = dict()
    grid_opts.setdefault("color", "black")
    grid_opts.setdefault("linestyle", "dotted")
    grid_opts.setdefault("linewidth", 1)
    ax.coords.grid(**grid_opts)

    return ax

def wcs_plot_world(ax, *args, **kwargs):
    """
    Wrapper for plotting lines/points in world coordinates.

    This mainly exists for convenience and not needing to remember how
    to use the transform argument.

    Arguments
    =========
    ax : WCSAxes instance
        The WCS axes, eg from wcs_subplot or wcsview

    Remaining positional and keyword arguments are passed ax.plot

    Returns
    =======
    the Line2D value(s) from ax.plot
    """
    kwargs["transform"] = ax.get_transform("world")
    return ax.plot(*args, **kwargs)

def wcsview(m, vmin=None, vmax=None, pmin=0.5, pmax=99.5, diverging=True,
        cmap=None, log=False, region="cmb", coord="C", projection="ZEA",
        cbar_extend='neither', cbar_ticks=None, cbar_minorticks=False, pad=.2, location='right',
        wcs=None, cbar=True, subplot=(), subplot_kwargs=None, **kwargs):
    """
    Plot a map using FITS WCS projections

    Arguments
    =========
    m : float, array-like or None
        An array containing the map,
        supports masked maps, see the `ma` function.
    vmin, vmax : float, optional
        The minimum and maximum range values. Override pmin, pmax.
        v-less forms for compatibility with healpy.
    pmin, pmax : float, optional
        Percentile minimum and maximum range values. Can also provide as a
        string with "%" to (v)min or (v)max.
    diverging : bool, optional
        Whether to set symmetric colour limits for diverging colour scale.
        Default based on region.
    cmap : string or colormap, optional
        The colormap to use. Default chosen based on diverging/not
    log : bool, optional
        Use a logarithmic rather than linear colour scale
    region : string, optional
        Determine projection parameters automatically based on the desired
        region to plot. Options include 'cmb' for main SPT3G region, or
        'latlon', 'rcw38', etc
    coord : sequence of character
        Either one of 'G', 'E' or 'C' to describe the coordinate
        system of the map
    projection : string, optional
        A projection name supported by basemap. Default: equal area azimuthal
        Must use the FITS WCS 3-letter name of the projection.
    wcs : WCS object, optional
        Specify an explicit WCS object. By default creates one based on
        region, coord, and projection.
        Overrides region, coord, and projection.
    subplot : tuple, optional
        Subplot specifier. Eg a tuple (nrows, ncols, index).
        Defaults to (1, 1, 1).  See pyplot.subplot
    subplot_kwargs : dict, optional
        Extra keyword arguments to pass to subplot (via wcs_subplot)

    Keyword arguments are passed to imshow.

    Returns
    =======
    (ax, im) : tuple of the axes (from wcs_subplot) and the image (from imshow)
    """
    # adjust arguments and default values
    if cmap is None:
        cmap = "merge" if diverging else "inferno"
        #cmap = plt.cm.RdBu_r if diverging else "inferno"
    bgcolor = "0.85"
    if isinstance(cmap, str):
        from matplotlib import cm
        cmap = copy.copy(cm.get_cmap(cmap))
        cmap.set_bad(bgcolor)
    # color scale
    if '%' in str(vmin):
        pmin = float(vmin.rstrip('%'))
    if '%' in str(vmax):
        pmax = float(vmax.rstrip('%'))
    if vmin is None or vmax is None:
        mask = mask_good(m, badinf=True)
        if vmin is None:
            vmin = np.percentile(m[mask], pmin)
        if vmax is None:
            vmax = np.percentile(m[mask], pmax)
    if diverging:
        val = np.max([np.abs(vmin), np.abs(vmax)])
        vmin = -val
        vmax= val
    if log:
        from matplotlib.colors import LogNorm
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    if subplot_kwargs is None:
        subplot_kwargs = dict()

    # create the WCS axes and plot
    #fig = plt.figure()
    if wcs is None:
        wcs = get_wcs_for_region(region=region, coord=coord,
                projection=projection)
    ax = wcs_subplot(wcs, subplot=subplot, **subplot_kwargs)
    mwcs = healpix_to_wcs(m, wcs)
    im = ax.imshow(mwcs, vmin=vmin, vmax=vmax, origin="lower",
            cmap=cmap, norm=norm, **kwargs)
    
    if cbar:
        figure = ax.get_figure()
        divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.1)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cb = figure.colorbar(im,ax=ax, location=location, label='$\mu$K', extend=cbar_extend)#, pad=pad, ax=ax, format='%.0e')
        #cb = figure.colorbar(im, cax=ax, label='$\mu$K')
        cax = figure.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        plt.colorbar(im, cax=cax)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
        if cbar_minorticks:
            cb.ax.minorticks_on()
        plt.sca(ax)     # reset current axis to main image
    else:
        cb = None

    im.set_clip_path(ax.coords.frame.patch)

    return ax, im

def make_cmap(data, name="new_cmap", N=256, register=False):
    """
    Make a linear segmented colour map from a table of data

    Arguments
    ---------
    data : shape (N, 3) array of RGB colours.
        The data for the cmap.
    name : sting
        Name for the cmap.
    N : int
        Number of RGB quantization levels.
    register : bool
        If True, will register the new cmap with matplotlib for get_cmap.

    Returns: the new colour map.
    """
    from matplotlib import colors, cm

    if np.max(data) > 1:
        data /= 256.
    cmap = colors.LinearSegmentedColormap.from_list(name, data, N=N)
    if register:
        cm.register_cmap(name=name, cmap=cmap)
    return cmap

def load_cmap(path, name=None, N=256, register=False, register_reverse=True):
    """
    Convenience function to load colour data and call make_cmap

    Arguments
    ---------
    path : string
        Path to the file. Should three columns of red, green, blue values.
    name : sting
        Name for the cmap. Will guess from filename.
    N : int
        Number of RGB quantization levels.
    register : bool
        If True, will register the new cmap with matplotlib for get_cmap.
    register_reverse : bool
        If this and register are True, will also register a cmap with
        colour order reversed and "_r" appended to name.

    Returns: the new colour map.
    """
    data = np.loadtxt(path)
    if not name:
        name = os.path.splitext(os.path.basename(path))[0]
        if name.endswith("_rgb"):
            name = name[:-4]
    cmap = make_cmap(data, name, register=register)
    if register and register_reverse:
        make_cmap(data[::-1,:], name+"_r", register=register)
    return cmap

def mapTQU2mapTEB(in_map):
    n_components = in_map.shape[0]
    nside = hp.npix2nside(in_map.shape[1])
    alms = hp.map2alm(in_map, pol=True)
    return np.array([in_map[0], hp.alm2map(alms[1], nside=nside), hp.alm2map(alms[2], nside=nside)])

def mask_bad(m, badval=hp.UNSEEN, rtol=1.e-5, atol=1.e-8,
             badnan=True, badinf=True):
    """Returns a bool array with ``True`` where m is close to badval,
    NaN or inf.

    Parameters
    ----------
    m : a map (may be a sequence of maps)
    badval : float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol : float, optional
        The relative tolerance
    atol : float, optional
        The absolute tolerance
    badnan : bool, optional
        If True, also mask NaN values
    badinf : bool, optional
        If True, also mask inf values

    Returns
    -------
    mask
      a bool array with the same shape as the input map, ``True`` where input map is
      close to badval, NaN or inf, and ``False`` elsewhere.

    See Also
    --------
    mask_good

    Examples
    --------
    >>> import healpy as hp
    >>> import numpy as np
    >>> m = np.arange(12.)
    >>> m[3] = hp.UNSEEN
    >>> m[4] = np.nan
    >>> mask_bad(m)
    array([False, False, False,  True,  True, False, False, False, False,
           False, False, False], dtype=bool)
    """
    m = np.asarray(m)
    mask = np.zeros_like(m, dtype=bool)
    if badnan:
        mask |= np.isnan(m)
    if badinf:
        mask |= np.isinf(m)
    mask[~mask] = hp.mask_bad(m[~mask], badval=badval, rtol=rtol, atol=atol)
    return mask

def smoothalm(alms, fwhm=0.0, sigma=None, beam=None, pol=True,
              mmax=None, verbose=True, inplace=True):
    """Smooth alm with a Gaussian symmetric beam or custom window function.

    Parameters
    ----------
    alms : array or sequence of 3 arrays
      Either an array representing one alm, or a sequence of arrays.
      See *pol* parameter.
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
      [in radians]
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.
      [in radians]
    beam : array or sequence of 3 arrays, optional
      If supplied, the beam function is applied instead of a Gaussian
      beam to each alm.
    pol : bool, optional
      If True, assumes input alms are TEB. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, apply spin 0 harmonic transform to each alm.
      (input can be any number of alms)
      If there is only one input alm, it has no effect. Default: True.
    mmax : None or int, optional
      The maximum m for alm. Default: mmax=lmax
    inplace : bool, optional
      If True, the alm's are modified inplace if they are contiguous arrays
      of type complex128. Otherwise, a copy of alm is made. Default: True.
    verbose : bool, optional
      If True prints diagnostic information. Default: True

    Returns
    -------
    alms : array or sequence of 3 arrays
      The smoothed alm. If alm[i] is a contiguous array of type complex128,
      and *inplace* is True the smoothing is applied inplace.
      Otherwise, a copy is made.
    """

    # make imports identical to healpy source for easy porting
    from healpy.sphtfunc import almxfl, Alm
    from healpy import cookbook as cb
    import numpy as np
    import six

    if beam is None:
        if sigma is None:
            sigma = fwhm / (2.*np.sqrt(2.*np.log(2.)))

        if verbose:
            print("Sigma is {0:f} arcmin ({1:f} rad) ".format(sigma*60*180/np.pi,sigma))
            print("-> fwhm is {0:f} arcmin".format(sigma*60*180/np.pi*(2.*np.sqrt(2.*np.log(2.)))))

    # Check alms
    if not cb.is_seq(alms):
        raise ValueError("alm must be a sequence")

    if sigma == 0 and beam is None:
        # nothing to be done
        return alms

    lonely = False
    if not cb.is_seq_of_seq(alms):
        alms = [alms]
        lonely = True

    # check beam
    if beam is not None:
        if not cb.is_seq(beam):
            raise ValueError("beam must be a sequence")
        if not lonely:
            if not cb.is_seq_of_seq(beam):
                beam = [beam]*len(alms)
            else:
                if len(beam) != len(alms):
                    raise ValueError("alm and beam shape mismatch")
        else:
            if cb.is_seq_of_seq(beam):
                raise ValueError("alm and beam shape mismatch")
            else:
                beam = [beam]

    # we have 3 alms -> apply smoothing to each map.
    # polarization has different B_l from temperature
    # exp{-[ell(ell+1) - s**2] * sigma**2/2}
    # with s the spin of spherical harmonics
    # s = 2 for pol, s=0 for temperature
    retalm = []
    for ialm, alm in enumerate(alms):
        lmax = Alm.getlmax(len(alm), mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given '
                            'mmax (len(alms[%d]) = %d).'%(ialm, len(alm)))
        if beam is None:
            ell = np.arange(lmax + 1.)
            s = 2 if ialm >= 1 and pol else 0
            fact = np.exp(-0.5 * (ell * (ell + 1) - s ** 2) * sigma ** 2)
        else:
            fact = beam[ialm]
        res = almxfl(alm, fact, mmax = mmax, inplace = inplace)
        retalm.append(res)
    # Test what to return (inplace/not inplace...)
    # Case 1: 1d input, return 1d output
    if lonely:
        return retalm[0]
    # case 2: 2d input, check if in-place smoothing for all alm's
    for i in six.moves.xrange(len(alms)):
        samearray = alms[i] is retalm[i]
        if not samearray:
            # Case 2a:
            # at least one of the alm could not be smoothed in place:
            # return the list of alm
            return retalm
    # Case 2b:
    # all smoothing have been performed in place:
    # return the input alms
    return alms


def make_Tonly_g3map(m, store=False, filename='Tonly_map.g3'):
    out_m = core.G3Frame(core.G3FrameType.Map)

    out_m['T'] = G3map['T'] 
    out_m['Q'] = 0*G3map['Q'] 
    out_m['U'] = 0*G3map['U'] 

    if store:
        core.G3Writer(filename)(out_m) 
        return out_m
    else:
        return out_m

def subtract_map_T2P_leakage(freq, G3map):
    freq = str(freq)
    freqs_dict = {'90':'100', '150':'143', '220':'217'}
    template = [fr for fr in core.G3File('/sptgrid/user/rgualtie/PlanckMockTonly/HFI_SkyMap_'+freqs_dict[freq]+'_2048_R3.01_fullmission_winter_large_field_Tonly/'+freq+'GHz/coadd_HFI_SkyMap_'+freqs_dict[freq]+'_2048_R3.01_fullmission_winter_large_field_Tonly_'+freq+'GHz_nstubs256.g3')][0]
    maps.RemoveWeights(template, zero_nans=True) 
    out_m = core.G3Frame(core.G3FrameType.Map)
    out_m['T'] = G3map['T'] 
    out_m['Q'] = G3map['Q'] - template['Q'] 
    out_m['U'] = G3map['U'] - template['U'] 

    return out_m

def map_deproject_T(G3map, freq, mask):
    #Find alphaQ and alphaU
    freq = str(freq)
    _, T = cs.spectrum_spice(G3map['T'], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=False, tolerance=1.e-7)
    _, crossQ = cs.spectrum_spice(G3map['T'], map2=G3map['Q'], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=False, tolerance=1.e-7)
    _, crossU = cs.spectrum_spice(G3map['T'], map2=G3map['U'], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=False, tolerance=1.e-7)
    noise = np.load('/sptgrid/user/rgualtie/spectra/'+freq+'GHz/noise/'+freq+'GHz_noise_spectra.npy')
    alphaQ = float(np.average(crossQ[2:12]/T[2:12], weights=1./noise[0, 2:12]))  #float(np.sum(crossQ[:21]/T[:21]))
    alphaU = float(np.average(crossU[2:12]/T[2:12], weights=1./noise[0, 2:12]))  #float(np.sum(crossU[:21]/T[:21]))
    print(alphaQ, alphaU)
    out_m = core.G3Frame(core.G3FrameType.Map)
    out_m['T'] = G3map['T'] 
    out_m['Q'] = G3map['Q'] - alphaQ * G3map['T'] 
    out_m['U'] = G3map['U'] - alphaU * G3map['T'] 
    out_m['Wpol'] = G3map['Wpol']
    return out_m, alphaQ, alphaU

def fit_foreground_template_min_spec(m, template, p0=[0.04], mask=None,
                                     joint_fit=False, smooth=False,
                                     fwhm=50.5, in_place=False,
                                     half_missions=True):
    '''
    Fits the template to the input map by minimizing the power in the 
    three TT (for I coeffient) and EE (for QU coefficients) bins 
    below l=100 (delta_l = 25, l_min=8)

    Right now, this only does the EE fit-- nothing with T is implemented.

    Args
    ----------
    m : float array
        The input map[s]. If multiple provided, will use all the cross
        spectra of those maps and exclude auto spectra. Should have 
        T->P leakage subtracted.
    template : list of float_array
        The foreground templates. Should be T->P leakage subtracted.

    Optional Args
    -------------
    p0 : list of float
        The initial guesses for the fits. Must be same length as template.
    mask : bool array
        The region to mask. Default None.
    smooth : bool
        Presmooth the map before fitting. Smooths both m and template. Default
        is False
    fwhm : float
        If smooth is True, then smooth with a beam with full-width half max
        fwhm in arcmin. Default is 50.5.
    in_place : bool
        If True, does not create a copy of the data to perform operations. This
        will alter the input maps.
    half_missions: bool
        If True, maps and templates are assumed to be half missions and are
        paired together in cross spectra. Len(maps) must equal len(templates).
        Only one template coefficient is fit.
    Returns
    -------
    fit_amp : float
        Returns the fit parameter. If joint_fit is False, returns two fit
        values, the first for I and the second for P. If the fit fails,
        returns np.nan
    '''
    m = np.array([np.asarray(m['T']), np.asarray(m['Q']), np.asarray(m['U'])])
    template = np.array([np.asarray(template['T']), np.asarray(template['Q']), np.asarray(template['U'])])
    if not in_place:
        m = np.copy(m)
        template = np.copy(template)
    
    if mask is None:
        mask = np.ones_like(template).astype(bool)
    # Need lists of templates and maps if only one provided
    if m.ndim == 2:
        m = [m] 
    if template.ndim == 2:
        template = [template]

    for i, m0 in enumerate(m):
        assert m0.shape[0] == 3, 'Maps do not all have IQU'
        for t0 in template:
            assert m0.shape == t0.shape, 'Template shape {} '.format(t0.shape)+\
                'differs from map shape {}'.format(m0.shape)
        m[i] *= mask
        if smooth:
            m[i] = hpt.smoothing(m0, fwhm=np.deg2rad(fwhm))

    for i, tp in enumerate(template):
        assert tp.shape[0] == 3, 'Templates do not all have IQU'
        template[i] *= mask
    
        if smooth:
            template[i] = hpt.smoothing(tp, fwhm=np.deg2rad(fwhm))

    if not half_missions:
        if len(p0) != len(template):
            print('Not enough alpha guesses compared to templates. '+
                  'Setting all the same.')
            p0 = np.repeat(p0[0], len(template))

    # Define the residual
    def residual(alpha, m, templates):
        print(alpha)
        dm = np.copy(m)
        if half_missions:
            for i, tp in enumerate(templates):
                dm[i] -= alpha[0]*tp
        else:    
            for i, tp in enumerate(templates):
                for j, m0 in enumerate(dm):
                    dm[j] -= alpha[i]*tp

        if len(dm) == 1:

            ells, dls = cs.spectrum_spice(dm[0], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=True, tolerance=1.e-7)

        else: # do all the crosses (excluding all autos)

            for n, comb in enumerate(combinations(dm, 2)):
                # running average
                if n == 0:
                    _, dls = cs.spectrum_spice(comb[0], comb[1], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=True, tolerance=1.e-7)
                else:
                    _, dls0 = cs.spectrum_spice(comb[0], comb[1], lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=True, tolerance=1.e-7)
                    dls += dls0 #This is corwin's modification 
        return(np.mean(1e6*dls[1][2:4]))
    
    # Solve for the fit parameter
    from scipy.optimize import minimize
    f_P = minimize(residual, p0, args=(m, template), method='Nelder-Mead')

    fit_amp = np.zeros((1, len(template))) * np.nan
    if f_P.success:
        fit_amp = f_P.x
    
    return fit_amp

def subtract_foregrounds(freq=None, G3map=None, mask=None):
    freq = str(freq)
    freqs_dict = {'90':'100', '150':'143', '220':'217'}
    template = [fr for fr in core.G3File('/sptgrid/user/rgualtie/PlanckMockForegrounds/p353-p100_dust_template.g3')][0]
#[fr for fr in core.G3File('/sptgrid/user/rgualtie/CommanderMock/HFI_CompMap_Foregrounds-commander_'+freqs_dict[freq]+'_2048_R3.00_winter_wide/'+freq+'GHz/coadd_HFI_CompMap_Foregrounds-commander_'+freqs_dict[freq]+'_2048_R3.00_winter_wide_'+freq+'GHz_nstubs256.g3')][0]
    #maps.RemoveWeights(template, zero_nans=True)
    alpha = fit_foreground_template_min_spec(G3map, template, p0=[0.], mask=mask, joint_fit=False, smooth=True, fwhm=.5, in_place=False)
    out_m = core.G3Frame(core.G3FrameType.Map)

    out_m['T'] = G3map['T'] - template['T'] * alpha[0]
    out_m['Q'] = G3map['Q'] - template['Q'] * alpha[0]
    out_m['U'] = G3map['U'] - template['U'] * alpha[0]
    out_m['Wpol'] = G3map['Wpol']
    return out_m, alpha

def rot_mat(alpha): #feed alpha in deg
    """
    Calculates the rotation matrix to null EB spectra
    Takes the rotation angle in degrees as an input
    The output is a numpy array 6x6
    """
    alpha = np.deg2rad(alpha)
    m = [[1, 0,0,0,0,0],
         [0,np.cos(2*alpha)**2, np.sin(2*alpha)**2,0,0,0],
         [0,np.sin(2*alpha)**2, np.cos(2*alpha)**2,0,0,0],
         [0,0,0,np.cos(2*alpha),0,0],
         [0,np.sin(4*alpha)/2, -np.sin(4*alpha)/2,0,1,0],
         [0,0,0,np.sin(2*alpha),0,1]]
    return np.array(m)

def rotate_spectra(dlsb, alphas, lmax=21):
    dlsb_r = []
    for b in range(lmax):
        dlsb_r.append(np.linalg.solve(rot_mat(alphas[b]), dlsb[:,b]))
    dlsb_r = np.array(dlsb_r).T
    #dlsb_r[:,0] = 0.
    return dlsb_r

def mixmat_cleaning(f, dlsb):
    f = str(f)
    mixmat = np.load('bincenter_mixmatrix_'+f+'GHz_TPmap.npy')
    return np.linalg.solve(mixmat, dlsb.flatten()).reshape(6,21) 

def load_spice_corr(filename='spice_cor.fits'):
    from astropy.io import fits
    corr = fits.open(filename)
    TT=np.zeros(6144)
    TQ=np.zeros(6144)
    TU=np.zeros(6144)
    for i in range(6144):
        TT[i] = corr[1].data[i][1]
        TQ[i] = corr[1].data[i][4]
        TU[i] = corr[1].data[i][5]
    return TT,TQ,TU

